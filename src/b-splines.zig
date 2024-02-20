const std = @import("std");

test "b-spline evaluation - 1st-degree, scalar output" {
    const CtrlPoint = f64;
    const Knot = f32;
    const context = makeFloatContext(CtrlPoint, Knot);

    const degree = 1;
    const knots = [_]Knot{ 0, 0, 3, 4, 4 };
    const ctrl_points = [knots.len - degree - 1]CtrlPoint{ 1, 2, 5 };

    const spline = BSpline1D(Knot){ .knots = &knots };

    // TODO add test case for duplicate interior knots

    // Result is zero outside of knots range
    try std.testing.expectEqual(null, spline.evalCtrlPointsAt(
        CtrlPoint,
        context,
        degree,
        &ctrl_points,
        -1,
    ));
    try std.testing.expectEqual(null, spline.evalCtrlPointsAt(
        CtrlPoint,
        context,
        degree,
        &ctrl_points,
        5,
    ));

    // Matches control points on knots
    try std.testing.expectApproxEqAbs(
        1,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            0,
        ) orelse unreachable,
        1e-5,
    );
    try std.testing.expectApproxEqAbs(
        2,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            3,
        ) orelse unreachable,
        1e-5,
    );
    try std.testing.expectApproxEqAbs(
        5,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            4,
        ) orelse unreachable,
        1e-5,
    );

    // Linearly interpolates between control points
    try std.testing.expectApproxEqAbs(
        1.25,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            0.75,
        ) orelse unreachable,
        1e-5,
    );
    try std.testing.expectApproxEqAbs(
        3.5,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            3.5,
        ) orelse unreachable,
        1e-5,
    );
}

test "b-spline evaluation - 2nd-degree, scalar output" {
    const CtrlPoint = f32;
    const Knot = f32;
    const context = makeFloatContext(CtrlPoint, Knot);

    const degree = 2;
    const knots = [_]Knot{ 0, 0, 0, 3, 4, 4, 4 };
    const ctrl_points = [knots.len - degree - 1]CtrlPoint{ 1, 2, 5, 1 };

    const spline = BSpline1D(Knot){ .knots = &knots };

    // TODO add test case for duplicate interior knots

    // Result is zero outside of knots range
    try std.testing.expectEqual(null, spline.evalCtrlPointsAt(
        CtrlPoint,
        context,
        degree,
        &ctrl_points,
        -1,
    ));
    try std.testing.expectEqual(null, spline.evalCtrlPointsAt(
        CtrlPoint,
        context,
        degree,
        &ctrl_points,
        5,
    ));

    // Matches control points on edge knots
    try std.testing.expectApproxEqAbs(
        1,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            0,
        ) orelse unreachable,
        1e-5,
    );
    try std.testing.expectApproxEqAbs(
        1,
        spline.evalCtrlPointsAt(
            CtrlPoint,
            context,
            degree,
            &ctrl_points,
            4,
        ) orelse unreachable,
        1e-5,
    );

    // Between knots
    {
        const x = 0.75;
        const k0 = 0.0;
        const k1 = 3.0;
        const k2 = 4.0;

        const N0_2 = 1;

        const N1_1 = ((k1 - x) / (k1 - k0)) * N0_2;
        const N1_2 = ((x - k0) / (k1 - k0)) * N0_2;

        const N2_0 = ((k1 - x) / (k1 - k0)) * N1_1;
        const N2_1 = ((((x - k0) / (k1 - k0)) * N1_1) + (((k2 - x) / (k2 - k0)) * N1_2));
        const N2_2 = ((x - k0) / (k2 - k0)) * N1_2;

        const exptResult = N2_0 * 1.0 + N2_1 * 2.0 + N2_2 * 5.0;
        try std.testing.expectEqual(@as(f32, 1.578125), exptResult);

        try std.testing.expectApproxEqAbs(
            @as(f32, exptResult),
            spline.evalCtrlPointsAt(
                CtrlPoint,
                context,
                degree,
                &ctrl_points,
                x,
            ) orelse unreachable,
            1e-7,
        );
    }
    {
        const x = 3.25;
        const k0 = 0.0;
        const k1 = 3.0;
        const k2 = 4.0;

        const N0_3 = 1;

        const N1_2 = ((k2 - x) / (k2 - k1)) * N0_3;
        const N1_3 = ((x - k1) / (k2 - k1)) * N0_3;

        const N2_1 = ((k2 - x) / (k2 - k0)) * N1_2;
        const N2_2 = ((((x - k0) / (k2 - k0)) * N1_2) + (((k2 - x) / (k2 - k1)) * N1_3));
        const N2_3 = ((x - k1) / (k2 - k1)) * N1_3;

        const exptResult = N2_1 * 2.0 + N2_2 * 5.0 + N2_3 * 1.0;
        try std.testing.expectEqual(@as(f32, 4.328125), exptResult);

        try std.testing.expectApproxEqAbs(
            @as(f32, exptResult),
            spline.evalCtrlPointsAt(
                CtrlPoint,
                context,
                degree,
                &ctrl_points,
                x,
            ) orelse unreachable,
            1e-7,
        );
    }
}

test "b-spline bases" {
    const T = f32;

    {
        const degree = 0;
        const knots = [6]T{ 0, 1, 3, 4, 4, 5 };
        const spline = BSpline1D(T){ .knots = &knots };

        try std.testing.expectEqual(null, spline.basesAt(degree, -1.00));
        {
            const result = spline.basesAt(degree, 0.25) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 0), result.index);
            try std.testing.expectEqualSlices(T, &[1]T{1}, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 2.00) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 1), result.index);
            try std.testing.expectEqualSlices(T, &[1]T{1}, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 3.75) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 2), result.index);
            try std.testing.expectEqualSlices(T, &[1]T{1}, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 4.50) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 4), result.index);
            try std.testing.expectEqualSlices(T, &[1]T{1}, &result.basis_values);
        }
        try std.testing.expectEqual(null, spline.basesAt(degree, 6.00));
    }

    {
        const degree = 1;
        const knots = [8]T{ 0, 0, 1, 3, 3, 3, 4, 5 };
        const spline = BSpline1D(T){ .knots = &knots };

        try std.testing.expectEqual(null, spline.basesAt(degree, -1.00));
        {
            const result = spline.basesAt(degree, 0.25) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 0), result.index);
            try std.testing.expectEqualSlices(T, &[2]T{ 0.75, 0.25 }, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 2.00) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 1), result.index);
            try std.testing.expectEqualSlices(T, &[2]T{ 0.50, 0.50 }, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 3.75) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 4), result.index);
            try std.testing.expectEqualSlices(T, &[2]T{ 0.25, 0.75 }, &result.basis_values);
        }
        {
            const result = spline.basesAt(degree, 4.50) orelse unreachable;
            try std.testing.expectEqual(@as(usize, 5), result.index);
            try std.testing.expectEqualSlices(T, &[2]T{ 0.50, 0.50 }, &result.basis_values);
        }
        try std.testing.expectEqual(null, spline.basesAt(degree, 6.00));
    }
}

// TODO add knot ops to context (enables use of Q64 for knots)
pub fn BSpline1D(comptime Knot: type) type {
    const CmpContext = void;
    const cmpContext: CmpContext = {};
    const lessThan = std.sort.asc(Knot);
    const cmp = struct {
        pub fn call(context: CmpContext, a: Knot, b: Knot) std.math.Order {
            if (lessThan(context, a, b)) {
                return std.math.Order.lt;
            } else if (a == b) {
                return std.math.Order.eq;
            } else {
                return std.math.Order.gt;
            }
        }
    }.call;

    return struct {
        const Self = @This();

        /// The knots of the spline. Must be sorted in ascending order.
        knots: []const Knot,

        // TODO add SIMD routine for evaluation of many values in bulk
        pub fn evalCtrlPointsAt(
            self: Self,
            comptime CtrlPoint: type,
            comptime context: anytype,
            comptime degree: comptime_int,
            ctrl_points: []const CtrlPoint,
            x: Knot,
        ) ?CtrlPoint {
            std.debug.assert(ctrl_points.len + degree + 1 == self.knots.len);

            const bases = self.basesAt(degree, x) orelse return null;

            var result: CtrlPoint = context.makeZero();
            for (bases.basis_values, ctrl_points[bases.index .. bases.index + degree + 1]) |b, cp| {
                result = context.add(result, context.mulKnot(cp, b));
            }
            return result;
        }

        pub fn basesAt(
            self: Self,
            comptime degree: comptime_int,
            x: Knot,
        ) ?struct { index: usize, basis_values: [degree + 1]Knot } {
            const x_i = self.knotIndexContaining(x) orelse return null;
            const index = std.math.sub(usize, x_i, degree + 1) catch return null;
            var basis_values: [degree + 1]Knot = undefined;
            self.basisSplinesAbout(degree, &basis_values, x_i, x);

            return .{ .index = index, .basis_values = basis_values };
        }

        fn basisSplinesAbout(
            self: Self,
            comptime degree: comptime_int,
            out: *[degree + 1]Knot,
            i_0: usize,
            x: Knot,
        ) void {
            var buffer: @TypeOf(out.*) = undefined;

            var buf_prev = out;
            var buf_next: @TypeOf(out) = &buffer;
            if (degree % 2 == 1) {
                std.mem.swap(@TypeOf(out), &buf_prev, &buf_next);
            }

            buf_prev[0] = 1;
            inline for (1..degree + 1) |d| {
                defer std.mem.swap(@TypeOf(out), &buf_prev, &buf_next);

                const bases_prev: *[d]Knot = @ptrCast(buf_prev);
                const bases_next: *[d + 1]Knot = @ptrCast(buf_next);

                bases_next[0] = 0;
                for (bases_prev, 0..) |basis_prev, j| {
                    const a = self.alpha(d, i_0 + j, x) orelse {
                        bases_next[j + 1] = 0;
                        continue;
                    };
                    bases_next[j] += basis_prev * (1 - a);
                    bases_next[j + 1] = basis_prev * a;
                }
            }

            return;
        }

        fn alpha(self: Self, comptime degree: comptime_int, index_last: usize, x: Knot) ?Knot {
            comptime std.debug.assert(degree > 0);
            // TODO change range checks into assertions, change calling
            // functions to avoid invalid sections of range
            // -> avoids unnecessary checks
            if (index_last >= self.knots.len or index_last < degree)
                return null;
            std.debug.assert(self.knots[index_last - degree] <= x);
            std.debug.assert(x <= self.knots[index_last]);

            const denominator = self.knots[index_last] - self.knots[index_last - degree];
            if (denominator == 0) return null;
            return (x - self.knots[index_last - degree]) / denominator;
        }

        fn knotIndexContaining(self: Self, x: Knot) ?usize {
            const i = bisect(Knot, x, self.knots, {}, cmp, .left);
            if (i == 0) {
                if (lessThan(cmpContext, x, self.knots[0])) {
                    // less than the smallest knot -> not contained by any range
                    return null;
                }
                // Make sure to use inner knots -> avoid any repeats in the front.
                return bisect(Knot, x, self.knots, {}, cmp, .right);
            }
            if (i >= self.knots.len) {
                // greater than the largest knot -> not contained by any range
                return null;
            }
            return i;
        }
    };
}

pub fn makeFloatContext(comptime CtrlPoint: type, comptime Knot: type) type {
    comptime switch (@typeInfo(CtrlPoint)) {
        .Float => {},
        else => std.debug.panic(
            "control point type must be a float type; instead got {}",
            .{@typeName(CtrlPoint)},
        ),
    };
    comptime switch (@typeInfo(Knot)) {
        .Float => {},
        else => std.debug.panic(
            "knot type must be a float type; instead got {}",
            .{@typeName(Knot)},
        ),
    };

    return struct {
        pub fn makeZero() CtrlPoint {
            return 0.0;
        }
        pub fn add(cp1: CtrlPoint, cp2: CtrlPoint) CtrlPoint {
            return cp1 + cp2;
        }
        pub fn mulKnot(ctrl_point: CtrlPoint, knot: Knot) CtrlPoint {
            return ctrl_point * @as(CtrlPoint, knot);
        }
    };
}

test bisect {
    const T = i32;
    const Context = void;
    const cmp = struct {
        pub fn call(context: Context, a: T, b: T) std.math.Order {
            if (std.sort.asc(T)(context, a, b)) {
                return std.math.Order.lt;
            } else if (a == b) {
                return std.math.Order.eq;
            } else {
                return std.math.Order.gt;
            }
        }
    }.call;
    const context: Context = {};

    const array = [_]T{ 1, 3, 5, 8, 9 };
    const x_existing: T = 3;
    const x_new: T = 7;

    try std.testing.expectEqual(bisect(T, x_existing, &array, context, cmp, .left), 1);
    try std.testing.expectEqual(bisect(T, x_existing, &array, context, cmp, .right), 2);
    try std.testing.expectEqual(bisect(T, x_new, &array, context, cmp, .left), 3);
    try std.testing.expectEqual(bisect(T, x_new, &array, context, cmp, .right), 3);
}

fn bisect(
    comptime T: type,
    key: anytype,
    items: []const T,
    context: anytype,
    comptime compareFn: fn (@TypeOf(context), @TypeOf(key), T) std.math.Order,
    comptime side: enum { left, right },
) usize {
    var left: usize = 0;
    var right: usize = items.len;

    while (left < right) {
        // Avoid overflowing in the midpoint calculation
        const mid = left + (right - left) / 2;

        switch (compareFn(context, key, items[mid])) {
            .lt => right = mid,
            .gt => left = mid + 1,
            .eq => switch (side) {
                .left => right = mid,
                .right => left = mid + 1,
            },
        }
    }

    return left;
}
