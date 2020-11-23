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
    coeff_x: Vec<[f32; 15]>,
    coeff_y: Vec<[[f32; 4]; 5]>,
}

fn horner(coeff: &[f32], x: f32) -> f32 {
    coeff.iter().fold(0., |acc, &c| acc * x + c)
}

fn indefinite_int(coeff_x: &[f32; 15], coeff_y: &[[f32; 4]; 5], x: f32, y: f32, dxdy: f32) -> f32 {
    assert!(x >= 0. && x <= 1. && y >= 0. && y <= 1.);
    y * (horner(&coeff_y[0], y) * horner(&coeff_x[0..5], x) +
        if dxdy == 0. { 0. } else {
            let ay = y * dxdy;
            let ay2 = ay * ay;
            let ay3 = ay2 * ay;
            let ay4 = ay2 * ay2;
            -horner(&coeff_y[1], y) * horner(&coeff_x[5..9], x) * ay
                + horner(&coeff_y[2], y) * horner(&coeff_x[9..12], x) * ay2
                - horner(&coeff_y[3], y) * horner(&coeff_x[12..14], x) * ay3
                + horner(&coeff_y[4], y) * horner(&coeff_x[14..15], x) * ay4
        })
}

fn x_coeffs(coeff_x: &[f32; 4]) -> [f32; 15] {
    [coeff_x[3] / 4., coeff_x[2] / 3., coeff_x[1] / 2., coeff_x[0], 0.,
        coeff_x[3], coeff_x[2], coeff_x[1], coeff_x[0],
        coeff_x[3] * 3., coeff_x[2] * 2., coeff_x[1],
        coeff_x[3] * 6., coeff_x[2] * 2.,
        coeff_x[3] * 6., ]
}

fn y_coeffs(coeff_y: &[f32; 4]) -> [[f32; 4]; 5] {
    [[coeff_y[3] / 4., coeff_y[2] / 3., coeff_y[1] / 2., coeff_y[0]],
        [coeff_y[3] / 20., coeff_y[2] / 12., coeff_y[1] / 6., coeff_y[0] / 2.],
        [coeff_y[3] / 120., coeff_y[2] / 60., coeff_y[1] / 24., coeff_y[0] / 6.],
        [coeff_y[3] / 840., coeff_y[2] / 360., coeff_y[1] / 120., coeff_y[0] / 24.],
        [coeff_y[3] / 6720., coeff_y[2] / 2520., coeff_y[1] / 720., coeff_y[0] / 120.], ]
}

// TODO: is there a faster way? (investigate whether approx recip is good enough)
fn recip(x: f32) -> f32 {
    x.recip()
}

impl Raster {
    pub fn new(w: usize, h: usize) -> Raster {
        let b = 0.;
        let c = 0.;
        let coeffs = [[0., 0., -c, b / 6. + c],
            [b / 6., b / 2. + c, -5. * b / 2. - 2. * c + 3., 3. * b / 2. + c - 2.],
            [1. - b / 3., 0., 2. * b + c - 3., -3. * b / 2. - c + 2.],
            [b / 6., -b / 2. - c, b / 2. + 2. * c, -b / 6. - c]];
        Raster {
            w: w,
            h: h,
            a: vec![0.0; w * h + 4],
            coeff_x: coeffs.iter().rev().map(x_coeffs).collect(),
            coeff_y: coeffs.iter().rev().map(y_coeffs).collect(),
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
            for j in 2..3 {
                for i in 2..3 {
                    let ci = &self.coeff_x[i];
                    let cj = &self.coeff_y[j];
                    let linestart_x0i = ((y + j) * self.w) as isize + x0i + i as isize;
                    if linestart_x0i < 0 {
                        continue; // oob index
                    }
                    if x1i <= x0i + 1 {
                        let top = indefinite_int(ci, cj, 1. - xnext.fract(), y1.fract(), -dxdy);
                        let btm = indefinite_int(ci, cj, 1. - x.fract(), y0.fract(), -dxdy);
                        let int = top - btm;
                        let full = indefinite_int(ci, cj, 1., y1.fract(), 0.) - indefinite_int(ci, cj, 1., y0.fract(), 0.);

                        self.a[linestart_x0i as usize] += dir * int;
                        self.a[linestart_x0i as usize + 1] += dir * (full - int);
                    } else {
                        let s = (x1 - x0).recip();

                        if x < xnext {
                            let mut y_inter = y0 + s * (1.0 - x.fract()) * dy;
                            let top = indefinite_int(ci, cj, 0., y_inter.fract(), -dxdy);
                            let btm = indefinite_int(ci, cj, 1. - x.fract(), y0.fract(), -dxdy);
                            let mut int = top - btm;
                            self.a[linestart_x0i as usize] += dir * int;
                            let mut rect_btm = indefinite_int(ci, cj, 1., y0.fract(), 0.);
                            for xi in x0i + 1..x1i - 1 {
                                let rect_top = indefinite_int(ci, cj, 1., y_inter.fract(), 0.);
                                let y_inter_new = y_inter + s * dy;
                                let top = indefinite_int(ci, cj, 0., y_inter_new.fract(), -dxdy);
                                let btm = indefinite_int(ci, cj, 1., y_inter.fract(), -dxdy);
                                let int_new = top - btm;
                                self.a[(linestart_x0i - x0i + xi) as usize] += dir * (rect_top - rect_btm + int_new - int);
                                rect_btm = rect_top;
                                int = int_new;
                                y_inter = y_inter_new;
                            }
                            let rect_top = indefinite_int(ci, cj, 1., y_inter.fract(), 0.);
                            let top = indefinite_int(ci, cj, 1. - xnext.fract(), y1.fract(), -dxdy);
                            let btm = indefinite_int(ci, cj, 1., y_inter.fract(), -dxdy);
                            let int_new = top - btm;
                            self.a[(linestart_x0i - x0i + x1i - 1) as usize] += dir * (rect_top - rect_btm + int_new - int);
                            let rect_top_final = indefinite_int(ci, cj, 1., y1.fract(), 0.);
                            self.a[(linestart_x0i - x0i + x1i) as usize] += dir * (rect_top_final - rect_top - int);
                        } else {
                            let mut y_inter = y1 - s * (1.0 - x.fract()) * dy;
                            let top = indefinite_int(ci, cj, 1. - x.fract(), y1.fract(), -dxdy);
                            let btm = indefinite_int(ci, cj, 0., y_inter.fract(), -dxdy);
                            let mut int = top - btm;
                            self.a[linestart_x0i as usize] += dir * int;
                            let mut rect_top = indefinite_int(ci, cj, 1., y1.fract(), 0.);
                            for xi in x0i + 1..x1i - 1 {
                                let rect_btm = indefinite_int(ci, cj, 1., y_inter.fract(), 0.);
                                let y_inter_new = y_inter - s * dy;
                                let top = indefinite_int(ci, cj, 1., y_inter.fract(), -dxdy);
                                let btm = indefinite_int(ci, cj, 0., y_inter_new.fract(), -dxdy);
                                let int_new = top - btm;
                                self.a[(linestart_x0i - x0i + xi) as usize] += dir * (rect_top - rect_btm + int_new - int);
                                rect_top = rect_btm;
                                int = int_new;
                                y_inter = y_inter_new;
                            }
                            let rect_btm = indefinite_int(ci, cj, 1., y_inter.fract(), 0.);
                            let top = indefinite_int(ci, cj, 1., y_inter.fract(), -dxdy);
                            let btm = indefinite_int(ci, cj, 1. - xnext.fract(), y0.fract(), -dxdy);
                            let int_new = top - btm;
                            self.a[(linestart_x0i - x0i + x1i - 1) as usize] += dir * (rect_top - rect_btm + int_new - int);
                            let rect_btm_final = indefinite_int(ci, cj, 1., y1.fract(), 0.);
                            self.a[(linestart_x0i - x0i + x1i) as usize] += dir * (rect_btm - rect_btm_final - int);
                        }
                    }
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

    pub fn dump(&self) -> String {
        return format!("coeffs: {:?} {:?}\ndata: {:?}", self.coeff_x, self.coeff_y, self.a.iter().filter(|&&x| x != 0.).collect::<Vec<_>>());
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
