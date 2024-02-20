const std = @import("std");

test "results in segfault" {
    try std.testing.expectEqual(null, Struct(f32).func(0));
}

pub fn Struct(comptime T: type) type {
    return struct {
        pub fn func(comptime n: comptime_int) ?struct { x: [n]T } {
            return null;
        }
    };
}
