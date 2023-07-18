function transfrom_magnet_to_spher (Blx, Bly, Blz)

Blz_new = double(Blx);
Bly_new = double(Blz);
Blx_new = double(Bly);
Bl_new = [Blx_new'; Bly_new'; Blz_new'];

z_new = double(x);
y_new = double(z);
x_new = double(y);

Bl_sph = zeros(length(Bl), 3);
for i = 1 : length(Bl)
    Bl_sph(i,:) = (cart2sphvec(Bl_new(:,i),...
        atan_angle(y_new(i),x_new(i))*180/pi,...
        atan_angle(z_new(i),sqrt(x_new(i)^2 + y_new(i)^2)*180/pi )))';
end
figure
plot(epoch, Bl_sph)
datetick;
legend('B_{az}','B_{el}','B_{r}')
grid on

end