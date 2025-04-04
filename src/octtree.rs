use std::{
    ops::Range,
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{body::Body, partition::Partition};
use ultraviolet::Vec3;

use rayon::prelude::*;

#[derive(Clone, Copy)]
pub struct Oct {
    pub center: Vec3,
    pub size: f32,
}

impl Oct {
    pub fn new_containing(bodies: &[Body]) -> Self {
        let mut min_x = f32::MAX;
        let mut min_y = f32::MAX;
        let mut max_x = f32::MIN;
        let mut max_y = f32::MIN;
        let mut min_z = f32::MAX;
        let mut max_z = f32::MIN;

        for body in bodies {
            min_x = min_x.min(body.pos.x);
            min_y = min_y.min(body.pos.y);
            max_x = max_x.max(body.pos.x);
            max_y = max_y.max(body.pos.y);
            min_z = min_z.min(body.pos.z);
            max_z = max_z.max(body.pos.z);
        }

        let center: Vec3 = Vec3::new(min_x + max_x, min_y + max_y, min_z + max_z) * 0.5;
        let size = ((max_x - min_x).max(max_y - min_y)).max(max_z - min_z);

        Self { center, size }
    }

    pub fn into_octant(mut self, octant: usize) -> Self {
        self.size *= 0.5;
        self.center.x += ((octant & 1) as f32 - 0.5) * self.size;
        self.center.y += (((octant >> 1) & 1) as f32 - 0.5) * self.size;
        self.center.z += (0.5 - (octant >> 2) as f32) * self.size;
        self
    }

    pub fn subdivide(&self) -> [Oct; 8] {
        [0, 1, 2, 3, 4, 5, 6, 7].map(|i| self.into_octant(i))
    }
}

#[derive(Clone)]
pub struct Node {
    pub children: usize,
    pub next: usize,
    pub pos: Vec3,
    pub mass: f32,
    pub oct: Oct,
    pub bodies: Range<usize>,
}

impl Node {
    pub const ZEROED: Self = Self {
        children: 0,
        next: 0,
        pos: Vec3 { x: 0.0, y: 0.0, z: 0.0 },
        mass: 0.0,
        oct: Oct {
            center: Vec3 { x: 0.0, y: 0.0, z: 0.0},
            size: 0.0,
        },
        bodies: 0..0,
    };

    pub fn new(next: usize, oct: Oct, bodies: Range<usize>) -> Self {        
        Self {    
            children: 0,
            next,
            pos: Vec3::zero(),
            mass: 0.0,
            oct,
            bodies,
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children == 0
    }

    pub fn is_branch(&self) -> bool {
        self.children != 0
    }

    pub fn is_empty(&self) -> bool {
        self.mass == 0.0
    }
}

pub struct Octtree {
    pub t_sq: f32,
    pub e_sq: f32,
    pub leaf_capacity: usize,
    pub thread_capacity: usize,
    pub atomic_len: AtomicUsize,
    pub nodes: Vec<Node>,
    pub parents: Vec<usize>,
}

impl Octtree {
    pub const ROOT: usize = 0;

    pub fn new(theta: f32, epsilon: f32, leaf_capacity: usize, thread_capacity: usize) -> Self {        
        Self {
            t_sq: theta * theta,
            e_sq: epsilon * epsilon,
            leaf_capacity,
            thread_capacity,
            atomic_len: 0.into(),
            nodes: Vec::new(),
            parents: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.atomic_len.store(0, Ordering::Relaxed);
    }

    pub fn subdivide(&mut self, node: usize, bodies: &mut [Body], range: Range<usize>) -> usize {
        let center = self.nodes[node].oct.center;
    
        let mut split = [range.start, 0, 0, 0, 0, 0, 0, 0, range.end];
    
        // Split bodies into 2 groups based on z-coordinate (above/below the center)
        let predicate_z = |body: &Body| body.pos.z < center.z;
        split[4] = split[0] + bodies[split[0]..split[8]].partition(predicate_z);
    
        // Split each z-partition group into 2 based on y-coordinate (above/below the center)
        let predicate_y = |body: &Body| body.pos.y < center.y;
        split[2] = split[0] + bodies[split[0]..split[4]].partition(predicate_y);
        split[6] = split[4] + bodies[split[4]..split[8]].partition(predicate_y);
    
        // Split each y-partition group into 2 based on x-coordinate (left/right of the center)
        let predicate_x = |body: &Body| body.pos.x < center.x;
        split[1] = split[0] + bodies[split[0]..split[2]].partition(predicate_x);
        split[3] = split[2] + bodies[split[2]..split[4]].partition(predicate_x);
        split[5] = split[4] + bodies[split[4]..split[6]].partition(predicate_x);
        split[7] = split[6] + bodies[split[6]..split[8]].partition(predicate_x);
    
        let len = self.atomic_len.fetch_add(1, Ordering::Relaxed);
        let children = len * 8 + 1;
        self.parents[len] = node;
        self.nodes[node].children = children;
    
        let nexts = [
            children + 1,
            children + 2,
            children + 3,
            children + 4,
            children + 5,
            children + 6,
            children + 7,
            self.nodes[node].next,
        ];
        let octs = self.nodes[node].oct.subdivide();
        for i in 0..8 {
            let bodies = split[i]..split[i + 1];
            self.nodes[children + i] = Node::new(nexts[i], octs[i], bodies);
        }
    
        children
    }

pub fn propagate(&mut self) {
    let len = self.atomic_len.load(Ordering::Relaxed);
    for &node in self.parents[..len].iter().rev() {
        let i = self.nodes[node].children;

        self.nodes[node].pos = self.nodes[i].pos
            + self.nodes[i + 1].pos
            + self.nodes[i + 2].pos
            + self.nodes[i + 3].pos
            + self.nodes[i + 4].pos
            + self.nodes[i + 5].pos
            + self.nodes[i + 6].pos
            + self.nodes[i + 7].pos;

        self.nodes[node].mass = self.nodes[i].mass
            + self.nodes[i + 1].mass
            + self.nodes[i + 2].mass
            + self.nodes[i + 3].mass
            + self.nodes[i + 4].mass
            + self.nodes[i + 5].mass
            + self.nodes[i + 6].mass
            + self.nodes[i + 7].mass;
    }

    self.nodes[0..len * 8 + 1].par_iter_mut().for_each(|node| {
        node.pos /= node.mass.max(f32::MIN_POSITIVE);
    });
}


    pub fn build(&mut self, bodies: &mut [Body]) {
        self.clear();

        let new_len = 8 * bodies.len() + 1024;
        self.nodes.resize(new_len, Node::ZEROED);
        self.parents.resize(new_len / 8, 0);

        let oct = Oct::new_containing(bodies);
        self.nodes[Self::ROOT] = Node::new(0, oct, 0..bodies.len());

        let (tx, rx) = crossbeam::channel::unbounded();
        tx.send(Self::ROOT).unwrap();

        let octtree_ptr = self as *mut Octtree as usize;
        let bodies_ptr = bodies.as_ptr() as usize;
        let bodies_len = bodies.len();

        let counter = AtomicUsize::new(0);
        rayon::broadcast(|_| {
            let mut stack = Vec::new();
            let octtree = unsafe { &mut *(octtree_ptr as *mut Octtree) };
            let bodies =
                unsafe { std::slice::from_raw_parts_mut(bodies_ptr as *mut Body, bodies_len) };

            while counter.load(Ordering::Relaxed) != bodies.len() {
                while let Ok(node) = rx.try_recv() {
                    let range = octtree.nodes[node].bodies.clone();
                    let len = octtree.nodes[node].bodies.len();

                    if range.len() >= octtree.thread_capacity {
                        let children = octtree.subdivide(node, bodies, range);
                        for i in 0..8 {
                            if !self.nodes[children + i].bodies.is_empty() {
                                tx.send(children + i).unwrap();
                            }
                        }
                        continue;
                    }

                    counter.fetch_add(len, Ordering::Relaxed);

                    stack.push(node);
                    while let Some(node) = stack.pop() {
                        let range = octtree.nodes[node].bodies.clone();
                        if range.len() <= octtree.leaf_capacity {
                            octtree.nodes[node].pos =
                                bodies[range.clone()].iter().map(|b| b.pos * b.mass).sum();
                            octtree.nodes[node].mass =
                                bodies[range.clone()].iter().map(|b| b.mass).sum();
                            continue;
                        }
                        let children = octtree.subdivide(node, bodies, range);
                        for i in 0..8 {
                            if !self.nodes[children + i].bodies.is_empty() {
                                stack.push(children + i);
                            }
                        }
                    }
                }
            }
        });

        self.propagate();
    }

    pub fn acc_pos(&self, pos: Vec3, bodies: &[Body]) -> Vec3 {
        let mut acc = Vec3::zero();

        let mut node = Self::ROOT;
        loop {
            let n = self.nodes[node].clone();

            let d = n.pos - pos;
            let d_sq = d.mag_sq();

            if n.oct.size * n.oct.size < d_sq * self.t_sq {
                let denom = (d_sq + self.e_sq) * d_sq.sqrt();
                acc += d * (n.mass / denom);

                if n.next == 0 {
                    break;
                }
                node = n.next;
            } else if n.is_leaf() {
                for i in n.bodies {
                    let body = &bodies[i];
                    let d = body.pos - pos;
                    let d_sq = d.mag_sq();

                    let denom = (d_sq + self.e_sq) * d_sq.sqrt();
                    acc += d * (body.mass / denom).min(f32::MAX);
                }

                if n.next == 0 {
                    break;
                }
                node = n.next;
            } else {
                node = n.children;
            }
        }

        acc
    }
    
    pub fn acc(&self, bodies: &mut Vec<Body>) {
        let bodies_ptr = std::ptr::addr_of_mut!(*bodies) as usize;

        bodies.par_iter_mut().for_each(|body| {
            let bodies = unsafe { &*(bodies_ptr as *const Vec<Body>) };
            body.acc = self.acc_pos(body.pos, bodies);
        });
    }
}