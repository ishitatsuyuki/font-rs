// Copyright 2015 Google Inc. All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//! An antialiased rasterizer for quadratic Beziers

use crate::accumulate::accumulate;
use crate::geom::Point;

// TODO: sort out crate structure. Right now we want this when compiling raster as a binary,
// but need it commented out when compiling showttf
//mod geom;

pub struct Raster {
    w: usize,
    h: usize,
    a: Vec<f32>,
}

fn base(x0: f32, x1: f32, y0: f32, y1: f32) -> f32 {
    assert!(0. <= x0 && x0 <= 1.);
    assert!(0. <= x1 && x1 <= 1.);
    assert!(-0.01 <= y0 && y0 <= 1.01, "y0: {}", y0);
    assert!(-0.01 <= y1 && y1 <= 1.01, "y1: {}", y1);
    assert!(y1 >= y0);
    (x0 * x0 * (3. * y0 + y1) + 2. * x0 * x1 * (y0 + y1) + x1 * x1 * (y0 + 3. * y1)) * (y1 - y0) / 24.
}

fn full(y1: f32, y2: f32) -> f32 {
    (y2 * y2 - y1 * y1) / 4.
}

// TODO: is there a faster way? (investigate whether approx recip is good enough)
fn recip(x: f32) -> f32 {
    x.recip()
}

impl Raster {
    pub fn new(w: usize, h: usize) -> Raster {
        Raster {
            w: w,
            h: h,
            a: vec![0.0; w * h + 4],
        }
    }

    pub fn draw_line(&mut self, p0: &Point, p1: &Point) {
        //println!("draw_line {} {}", p0, p1);
        if (p0.y - p1.y).abs() <= f32::EPSILON {
            return;
        }
        let (dir, p0, p1) = if p0.y < p1.y {
            (1.0, p0, p1)
        } else {
            (-1.0, p1, p0)
        };
        let dxdy = (p1.x - p0.x) / (p1.y - p0.y);
        let dydx = (p1.y - p0.y) / (p1.x - p0.x);
        let mut x = p0.x;
        let y0 = p0.y as usize; // note: implicit max of 0 because usize (TODO: really true?)
        if p0.y < 0.0 {
            x -= p0.y * dxdy;
        }
        for y in y0..self.h.min(p1.y.ceil() as usize) {
            let y0 = (y as f32).max(p0.y);
            let y1 = ((y + 1) as f32).min(p1.y);
            let dy = y1 - y0;
            let xnext = x + dxdy * dy;
            let (x0, x1) = if x < xnext { (x, xnext) } else { (xnext, x) };
            let x0floor = x0.floor();
            let x0i = x0floor as isize;
            let x1ceil = x1.ceil();
            let x1i = x1ceil as isize;
            let linestart_x0i = (y * self.w) as isize + x0i;
            let linestart_x0i_1 = linestart_x0i + self.w as isize;
            if linestart_x0i < 0 {
                continue; // oob index
            }
            if x1i <= x0i + 1 {
                let x0f = x - x0floor;
                let x1f = xnext - x0floor;
                let y0f = y0.fract();
                let y1f = y1 - y0.floor();
                let p00 = base(x0f, x1f, y0f, y1f);
                let p10 = base(1. - x0f, 1. - x1f, y0f, y1f);
                let p01 = base(x1f, x0f, 1. - y1f, 1. - y0f);
                let p11 = base(1. - x1f, 1. - x0f, 1. - y1f, 1. - y0f);
                let full0 = full(y0f, y1f);
                let full1 = full(1. - y1f, 1. - y0f);
                self.a[linestart_x0i_1 as usize] += dir * p10;
                self.a[linestart_x0i_1 as usize + 1] += dir * (2. * full0 - p00 - p10);
                self.a[linestart_x0i_1 as usize + 2] += dir * p00;
                self.a[linestart_x0i as usize] += dir * p11;
                self.a[linestart_x0i as usize + 1] += dir * (2. * full1 - p01 - p11);
                self.a[linestart_x0i as usize + 2] += dir * p01;
            } else {
                for xi in x0i..x1i {
                    // TODO: not left and right
                    let linestart_x0i = (y * self.w) as isize + xi;
                    let linestart_x0i_1 = linestart_x0i + self.w as isize;
                    let xl = x.max(xi as f32).min((xi + 1) as f32);
                    let xr = xnext.max(xi as f32).min((xi + 1) as f32);
                    let x0f = xl - xi as f32;
                    let x1f = xr - xi as f32;
                    let y0f = y0.fract() + (xl - x) * dydx;
                    let y1f = y0.fract() + (xr - x) * dydx;
                    let p00 = base(x0f, x1f, y0f, y1f);
                    let p10 = base(1. - x0f, 1. - x1f, y0f, y1f);
                    let p01 = base(x1f, x0f, 1. - y1f, 1. - y0f);
                    let p11 = base(1. - x1f, 1. - x0f, 1. - y1f, 1. - y0f);
                    let full0 = full(y0f, y1f);
                    let full1 = full(1. - y1f, 1. - y0f);
                    self.a[linestart_x0i_1 as usize] += dir * p10;
                    self.a[linestart_x0i_1 as usize + 1] += dir * (2. * full0 - p00 - p10);
                    self.a[linestart_x0i_1 as usize + 2] += dir * p00;
                    self.a[linestart_x0i as usize] += dir * p11;
                    self.a[linestart_x0i as usize + 1] += dir * (2. * full1 - p01 - p11);
                    self.a[linestart_x0i as usize + 2] += dir * p01;
                }
            }
            x = xnext;
        }
    }

    pub fn draw_quad(&mut self, p0: &Point, p1: &Point, p2: &Point) {
        //println!("draw_quad {} {} {}", p0, p1, p2);
        let devx = p0.x - 2.0 * p1.x + p2.x;
        let devy = p0.y - 2.0 * p1.y + p2.y;
        let devsq = devx * devx + devy * devy;
        if devsq < 0.333 {
            self.draw_line(p0, p2);
            return;
        }
        let tol = 3.0;
        let n = 1 + (tol * (devx * devx + devy * devy)).sqrt().sqrt().floor() as usize;
        //println!("n = {}", n);
        let mut p = *p0;
        let nrecip = recip(n as f32);
        let mut t = 0.0;
        for _i in 0..n - 1 {
            t += nrecip;
            let pn = Point::lerp(t, &Point::lerp(t, p0, p1), &Point::lerp(t, p1, p2));
            self.draw_line(&p, &pn);
            p = pn;
        }
        self.draw_line(&p, p2);
    }


    /*
    fn get_bitmap_fancy(&self) -> Vec<u8> {
        let mut acc = 0.0;
        // This would translate really well to SIMD
        self.a[0..self.w * self.h].iter().map(|&a| {
            acc += a;
            (255.0 * acc.abs().min(1.0)) as u8
            //(255.5 * (0.5 + 0.4 * acc)) as u8
        }).collect()
    }
*/

    pub fn get_bitmap(&self) -> Vec<u8> {
        accumulate(&self.a[0..self.w * self.h])
    }
}
