use std::array::IntoIter;

use freetype::freetype::*;
use kurbo::PathEl;
use kurbo::Point;

use crate::raster::Raster;

struct Context {
    raster: Raster,
    last_point: Option<Point>,
    height: i64,
}

unsafe fn vec_to_point(ctx: &Context, v: *const FT_Vector) -> Point {
    // TODO: better logic to handle padding
    Point::new((*v).x as f64 / 64., (ctx.height - 64 - (*v).y) as f64 / 64.)
}

unsafe extern "C" fn move_to(to: *const FT_Vector,
                             user: *mut ::std::os::raw::c_void) -> i32 {
    let ctx = &mut *(user as *mut Context);
    ctx.last_point = Some(vec_to_point(ctx, to));
    0
}

unsafe extern "C" fn line_to(to: *const FT_Vector,
                             user: *mut ::std::os::raw::c_void) -> i32 {
    let ctx = &mut *(user as *mut Context);
    let point = vec_to_point(ctx, to);
    if let Some(last) = ctx.last_point {
        ctx.raster.draw_path(IntoIter::new([PathEl::MoveTo(last), PathEl::LineTo(point)]))
    }
    ctx.last_point = Some(point);
    0
}

unsafe extern "C" fn conic_to(control1: *const FT_Vector, to: *const FT_Vector,
                              user: *mut ::std::os::raw::c_void) -> i32 {
    let ctx = &mut *(user as *mut Context);
    let control1 = vec_to_point(ctx, control1);
    let to = vec_to_point(ctx, to);
    if let Some(last) = ctx.last_point {
        ctx.raster.draw_path(IntoIter::new([PathEl::MoveTo(last), PathEl::QuadTo(control1, to)]))
    }
    ctx.last_point = Some(to);
    0
}

unsafe extern "C" fn cubic_to(control1: *const FT_Vector, control2: *const FT_Vector, to: *const FT_Vector,
                              user: *mut ::std::os::raw::c_void) -> i32 {
    let ctx = &mut *(user as *mut Context);
    let control1 = vec_to_point(ctx, control1);
    let control2 = vec_to_point(ctx, control2);
    let to = vec_to_point(ctx, to);
    if let Some(last) = ctx.last_point {
        ctx.raster.draw_path(IntoIter::new([PathEl::MoveTo(last), PathEl::CurveTo(control1, control2, to)]))
    }
    ctx.last_point = Some(to);
    0
}


pub fn raster_render(width: usize, height: usize, outline: &mut FT_Outline) -> Vec<u8> {
    let raster = Raster::new(width, height);
    let mut context = Context {
        raster,
        last_point: None,
        height: (height << 6) as _,
    };
    let funcs = FT_Outline_Funcs_ {
        move_to: Some(move_to),
        line_to: Some(line_to),
        conic_to: Some(conic_to),
        cubic_to: Some(cubic_to),
        shift: 0,
        delta: 0,
    };
    unsafe {
        FT_Outline_Decompose(outline, &funcs, &mut context as *mut _ as *mut _);
    }
    context.raster.get_bitmap()
}
