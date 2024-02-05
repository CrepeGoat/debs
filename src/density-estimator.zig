const std = @import("std");
const bspline = @import("b-splines.zig");

pub fn PdfEstimateBuilder(comptime T: type) type {
    const degree = 3;
    return struct {
        const Self = @This();
        const Knot = f64;

        samples: []const T,
        knots: std.ArrayList(Knot),
        spline: bspline.BSpline1D(Knot),
        spline_coeffs: std.ArrayList(T),
        spline_coeffs_prev: std.ArrayList(T),

        pub fn init(
            allocator: std.mem.Allocator,
            samples: []const T,
        ) std.Allocator.Error!Self {
            const coeffs_count = @min(samples.len, 128);
            const knots_count = coeffs_count + degree + 1;
            // const fudge_factor = 0.01;

            const knots = try std.ArrayList(Knot).initCapacity(allocator, knots_count);
            const spline = bspline.BSpline1D(Knot).initUniformIn(degree, &knots.items);
            var coeffs = try std.ArrayList(T).initCapacity(allocator, coeffs_count);
            coeffs.appendNTimesAssumeCapacity(0, coeffs_count);
            const coeffs_prev = try std.ArrayList(T).initCapacity(allocator, coeffs_count);

            return .{
                .samples = samples,
                .knots = knots,
                .spline = spline,
                .spline_coeffs = coeffs,
                .spline_coeffs_prev = coeffs_prev,
            };
        }

        pub fn deinit(self: Self) void {
            self.knots.deinit();
            self.spline_coeffs.deinit();
            self.spline_coeffs_prev.deinit();
        }

        pub fn incrementCoeffs(self: *Self) void {
            _ = self; // autofix

        }

        fn evalError(self: Self) void {
            for (self.samples) |sample| {
                const bases = self.spline.basesAt(degree, sample);
                _ = bases; // autofix
            }
        }

        fn estimateAt(self: Self, x: Knot) T {
            return self.spline.evalCtrlPointsAt(T, {}, degree, self.coeffs, x);
        }

        fn makeInitCoeffs(coeffs: []T) void {
            _ = coeffs; // autofix

        }
    };
}
