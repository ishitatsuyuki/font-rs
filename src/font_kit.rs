use ::font_kit::outline::OutlineSink;

struct OutlineRasterizer {

}

impl OutlineSink for OutlineRasterizer {
    fn move_to(&mut self, to: Vector2F) {
        unimplemented!()
    }

    fn line_to(&mut self, to: Vector2F) {
        unimplemented!()
    }

    fn quadratic_curve_to(&mut self, ctrl: Vector2F, to: Vector2F) {
        unimplemented!()
    }

    fn cubic_curve_to(&mut self, ctrl: LineSegment2F, to: Vector2F) {
        unimplemented!()
    }

    fn close(&mut self) {
        unimplemented!()
    }
}