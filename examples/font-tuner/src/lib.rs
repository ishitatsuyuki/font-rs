use std::iter;

use wasm_bindgen::prelude::*;

use font_rs::font::parse;
use font_rs::geom::Affine;
use font_rs::raster::Raster;

mod utils;

// A macro to provide `println!(..)`-style syntax for `console.log` logging.
macro_rules! log {
    ( $( $t:tt )* ) => {
        web_sys::console::log_1(&format!( $( $t )* ).into());
    }
}

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

fn trc(gamma: f32, luminance: f32) -> f32 {
    match gamma {
        _ if gamma == 0.0 => {
            //The magic numbers are derived from the sRGB specification.
            //See http://www.color.org/chardata/rgb/srgb.xalter .
            if luminance <= 0.04045 {
                luminance / 12.92
            } else {
                ((luminance + 0.055) / 1.055).powf(2.4)
            }
        }
        _ => luminance.powf(gamma),
    }
}

fn inverse_trc(gamma: f32, luma: f32) -> f32 {
    match gamma {
        _ if gamma == 0.0 => {
            //The magic numbers are derived from the sRGB specification.
            //See http://www.color.org/chardata/rgb/srgb.xalter .
            if luma <= 0.0031308 {
                luma * 12.92
            } else {
                1.055 * luma.powf(1. / 2.4) - 0.055
            }
        }
        _ => luma.powf(gamma.recip()),
    }
}

#[wasm_bindgen]
pub fn render(font_bytes: &[u8], text: &str) -> Box<[u8]> {
    utils::set_panic_hook();
    let gamma = 0.0;
    let gamma_lut: Vec<_> = (0..256).map(|x| x as f32 / 255.).map(|x| (inverse_trc(gamma, x) * 255. + 0.5) as u8).collect();
    let size = 12;
    let font = parse(font_bytes).unwrap();
    let mut raster = Raster::new(512, 512);
    let v_metrics = font.get_v_metrics(size).unwrap();
    let scale = font.scale(size);
    let mut pen = Affine::new(scale, 0., 0., -scale, 0.5, 0.5 + v_metrics.ascent);
    for c in text.chars() {
        if c == '\n' {
            pen.f += (v_metrics.ascent - v_metrics.descent + v_metrics.line_gap);
            pen.e = 0.5;
            continue;
        }
        let glyph = font.lookup_glyph_id(c as u32).unwrap();
        let h_metrics = font.get_h_metrics(glyph, size).unwrap();
        let glyph = font.get_glyph(glyph).unwrap();
        font.render_glyph_inner(&mut raster, &pen, &glyph);
        pen.e += h_metrics.advance_width;
    }
    let bitmap = raster.get_bitmap();
    bitmap.into_iter().flat_map(|x| iter::repeat(gamma_lut[x as usize]).take(3).chain(iter::once(255))).collect::<Vec<_>>().into_boxed_slice()
}
