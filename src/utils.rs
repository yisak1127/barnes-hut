use crate::body::Body;
use ultraviolet::Vec3;

pub fn uniform_disc(n: usize) -> Vec<Body> {
    fastrand::seed(0);
    let inner_radius = 25.0;
    let outer_radius = (n as f32).sqrt() * 5.0;

    let mut bodies: Vec<Body> = Vec::with_capacity(n);

    let m = 1e6;
    let center = Body::new(Vec3::zero(), Vec3::zero(), m as f32, inner_radius);
    bodies.push(center);

    while bodies.len() < n {
        let a = fastrand::f32() * std::f32::consts::TAU;
        let b = fastrand::f32() * std::f32::consts::TAU/2.0;
        let (sin, cos) = a.sin_cos();
        let (sin_phi, cos_phi) = b.sin_cos();
        let t = inner_radius / outer_radius;
        let r = fastrand::f32() * (1.0 - t * t) + t * t;
        let pos = Vec3::new(cos*sin_phi, sin*sin_phi, cos_phi) * outer_radius * r.sqrt();
        let vel = Vec3::new(sin*sin_phi, -cos*sin_phi, sin_phi);
        let mass = 1.0f32;
        let radius = mass.cbrt();

        bodies.push(Body::new(pos, vel, mass, radius));
    }

    bodies.sort_by(|a, b| a.pos.mag_sq().total_cmp(&b.pos.mag_sq()));
    let mut mass = 0.0;
    for i in 0..n {
        mass += bodies[i].mass;
        if bodies[i].pos == Vec3::zero() {
            continue;
        }

        let v = (mass / bodies[i].pos.mag()).sqrt();
        bodies[i].vel *= v;
    }

    bodies
}
