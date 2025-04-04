use std::{
    f32::consts::{PI, TAU},
    sync::atomic::{AtomicBool, Ordering},
};

use crate::{
    body::Body,
    octtree::{Node, Octtree},
};

use quarkstrom::{egui, winit::event::VirtualKeyCode, winit_input_helper::WinitInputHelper};

use palette::{rgb::Rgba, Hsluv, IntoColor};
use ultraviolet::Vec2;
use ultraviolet::Vec3;
use ultraviolet::Mat3;

use once_cell::sync::Lazy;
use parking_lot::Mutex;

pub static PAUSED: Lazy<AtomicBool> = Lazy::new(|| false.into());
pub static UPDATE_LOCK: Lazy<Mutex<bool>> = Lazy::new(|| Mutex::new(false));

pub static BODIES: Lazy<Mutex<Vec<Body>>> = Lazy::new(|| Mutex::new(Vec::new()));
pub static OCTTREE: Lazy<Mutex<Vec<Node>>> = Lazy::new(|| Mutex::new(Vec::new()));

pub static SPAWN: Lazy<Mutex<Vec<Body>>> = Lazy::new(|| Mutex::new(Vec::new()));

pub struct Renderer {
    pos: Vec3,
    scale: f32,

    rotation_x: f32, // Rotation around the X-axis (pitch)
    rotation_y: f32, // Rotation around the Y-axis (yaw)

    settings_window_open: bool,

    show_bodies: bool,
    show_octtree: bool,

    depth_range: (usize, usize),

    spawn_body: Option<Body>,
    angle: Option<f32>,
    total: Option<f32>,

    confirmed_bodies: Option<Body>,

    bodies: Vec<Body>,
    octtree: Vec<Node>,
}

impl quarkstrom::Renderer for Renderer {
    fn new() -> Self {
        Self {
            pos: Vec3::zero(),
            scale: 3600.0,

            rotation_x: 0.0, // Initialize rotations to 0
            rotation_y: 0.0,

            settings_window_open: false,

            show_bodies: true,
            show_octtree: false,

            depth_range: (0, 0),

            spawn_body: None,
            angle: None,
            total: None,

            confirmed_bodies: None,

            bodies: Vec::new(),
            octtree: Vec::new(),
        }
    }

    fn input(&mut self, input: &WinitInputHelper, width: u16, height: u16) {
        self.settings_window_open ^= input.key_pressed(VirtualKeyCode::E);

        if input.key_pressed(VirtualKeyCode::Space) {
            let val = PAUSED.load(Ordering::Relaxed);
            PAUSED.store(!val, Ordering::Relaxed)
        }

        if let Some((mx, my)) = input.mouse() {
            // Scroll steps to double/halve the scale
            let steps = 5.0;

            // Modify input
            let zoom = (-input.scroll_diff() / steps).exp2();

            // Screen space -> view space
            let target =
                Vec2::new(mx * 2.0 - width as f32, height as f32 - my * 2.0) / height as f32;

            // Move view position based on target
            self.pos += (target * self.scale * (1.0 - zoom)).into();

            // Zoom
            self.scale *= zoom;
        }

        if input.mouse_held(0) {
            let (mdx, mdy) = input.mouse_diff();
            self.rotation_y += mdx / height as f32 * self.scale * 0.01;
            self.rotation_x += mdy / height as f32 * self.scale * 0.01;
        }

        // Grab
        if input.mouse_held(2) {
            let (mdx, mdy) = input.mouse_diff();
            self.pos.x -= mdx / height as f32 * self.scale * 2.0;
            self.pos.y += mdy / height as f32 * self.scale * 2.0;
        }

        let world_mouse = || -> Vec3 {
            let (mx, my) = input.mouse().unwrap_or_default();
            let mut mouse = Vec2::new(mx, my);
            mouse *= 2.0 / height as f32;
            mouse.y -= 1.0;
            mouse.y *= -1.0;
            mouse.x -= width as f32 / height as f32;
            (mouse * self.scale + self.pos.into()).into()
        };

        if input.mouse_pressed(1) {
            let mouse = world_mouse();
            self.spawn_body = Some(Body::new(mouse, Vec3::zero(), 1.0, 1.0));
            self.angle = None;
            self.total = Some(0.0);
        } else if input.mouse_held(1) {
            if let Some(body) = &mut self.spawn_body {
                let mouse = world_mouse();
                if let Some(angle) = self.angle {
                    let d = mouse - body.pos;
                    let angle2 = d.y.atan2(d.x);
                    let a = angle2 - angle;
                    let a = (a + PI).rem_euclid(TAU) - PI;
                    let total = self.total.unwrap() - a;
                    body.mass = (total / TAU).exp2();
                    self.angle = Some(angle2);
                    self.total = Some(total);
                } else {
                    let d = mouse - body.pos;
                    let angle = d.y.atan2(d.x);
                    self.angle = Some(angle);
                }
                body.radius = body.mass.cbrt();
                body.vel = mouse - body.pos;
            }
        } else if input.mouse_released(1) {
            self.confirmed_bodies = self.spawn_body.take();
        }
    }

    fn render(&mut self, ctx: &mut quarkstrom::RenderContext) {
        {
            let mut lock = UPDATE_LOCK.lock();
            if *lock {
                std::mem::swap(&mut self.bodies, &mut BODIES.lock());
                std::mem::swap(&mut self.octtree, &mut OCTTREE.lock());
            }
            if let Some(body) = self.confirmed_bodies.take() {
                self.bodies.push(body);
                SPAWN.lock().push(body);
            }
            *lock = false;
        }
        let rotation_matrix_x = Mat3::from_rotation_x(self.rotation_x);
        let rotation_matrix_y = Mat3::from_rotation_y(self.rotation_y);
        let rotation_matrix = rotation_matrix_y * rotation_matrix_x;

        ctx.clear_circles();
        ctx.clear_lines();
        ctx.clear_rects();
        ctx.set_view_pos(self.pos.into());
        ctx.set_view_scale(self.scale);

        if !self.bodies.is_empty() {
            if self.show_bodies {
                for i in 0..self.bodies.len() {
                    let rotated_pos = rotation_matrix * self.bodies[i].pos;
                    ctx.draw_circle(rotated_pos.into(), self.bodies[i].radius, [0xff; 4]);
                }
            }

            if let Some(body) = &self.confirmed_bodies {
                ctx.draw_circle(body.pos.into(), body.radius, [0xff; 4]);
                ctx.draw_line(body.pos.into(), (body.pos + body.vel).into(), [0xff; 4]);
            }

            if let Some(body) = &self.spawn_body {
                ctx.draw_circle(body.pos.into(), body.radius, [0xff; 4]);
                ctx.draw_line(body.pos.into(), (body.pos + body.vel).into(), [0xff; 4]);
            }
        }

        if self.show_octtree && !self.octtree.is_empty() {
            let mut depth_range = self.depth_range;
            if depth_range.0 >= depth_range.1 {
                let mut stack = Vec::new();
                stack.push((Octtree::ROOT, 0));

                let mut min_depth = usize::MAX;
                let mut max_depth = 0;
                while let Some((node, depth)) = stack.pop() {
                    let node = &self.octtree[node];

                    if node.is_leaf() {
                        if depth < min_depth {
                            min_depth = depth;
                        }
                        if depth > max_depth {
                            max_depth = depth;
                        }
                    } else {
                        for i in 0..4 {
                            stack.push((node.children + i, depth + 1));
                        }
                    }
                }

                depth_range = (min_depth, max_depth);
            }
            let (min_depth, max_depth) = depth_range;

            let mut stack = Vec::new();
            stack.push((Octtree::ROOT, 0));
            while let Some((node, depth)) = stack.pop() {
                let node = &self.octtree[node];

                if node.is_branch() && depth < max_depth {
                    for i in 0..4 {
                        stack.push((node.children + i, depth + 1));
                    }
                } else if depth >= min_depth {
                    let oct = node.oct;
                    let half = Vec3::new(0.5, 0.5, 0.5) * oct.size;
                    let min = oct.center - half;
                    let max = oct.center + half;

                    let t = ((depth - min_depth + !node.is_empty() as usize) as f32)
                        / (max_depth - min_depth + 1) as f32;

                    let start_h = -100.0;
                    let end_h = 80.0;
                    let h = start_h + (end_h - start_h) * t;
                    let s = 100.0;
                    let l = t * 100.0;

                    let c = Hsluv::new(h, s, l);
                    let rgba: Rgba = c.into_color();
                    let color = rgba.into_format().into();

                    ctx.draw_rect(min.into(), max.into(), color);
                }
            }
        }
    }

    fn gui(&mut self, ctx: &quarkstrom::egui::Context) {
        egui::Window::new("")
            .open(&mut self.settings_window_open)
            .show(ctx, |ui| {
                ui.checkbox(&mut self.show_bodies, "Show Bodies");
                ui.checkbox(&mut self.show_octtree, "Show octtree");
                if self.show_octtree {
                    let range = &mut self.depth_range;
                    ui.horizontal(|ui| {
                        ui.label("Depth Range:");
                        ui.add(egui::DragValue::new(&mut range.0).speed(0.05));
                        ui.label("to");
                        ui.add(egui::DragValue::new(&mut range.1).speed(0.05));
                    });
                }
            });
    }
}
