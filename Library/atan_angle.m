% returns angle in radians
function angle = atan_angle (opposite_side,adjusent_side)
angle = atan(opposite_side./adjusent_side);
angle (opposite_side >= 0  & adjusent_side < 0) = pi + angle (opposite_side >= 0  & adjusent_side < 0);
angle (opposite_side < 0  & adjusent_side <= 0)= pi + angle (opposite_side < 0  & adjusent_side <= 0);
angle (opposite_side < 0  & adjusent_side > 0)= 2*pi + angle (opposite_side < 0  & adjusent_side > 0);