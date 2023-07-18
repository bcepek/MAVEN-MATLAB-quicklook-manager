function [SW_time, SW_Bx, SW_By, SW_Bz, SW_v, SW_n_p, SW_dyn_press] = read_Halekas_SW 

UT_num = 1;
n_p_num = 2;
n_He_num = 3;
v_num = 4;
v_x_num = 5;
v_y_num = 6;
v_z_num = 7;
T_num = 8;
Bx_num = 9;
By_num = 10;
Bz_num = 11;
% save ('\\193.232.6.100\общие\Владимир\MAVEN upstream parameters from Halekas 01-08-2017\drivers_merge_l2.mat', 'Halekas_SW')
load('\\193.232.6.100\общие\Владимир\MAVEN upstream parameters from Halekas 01-08-2017\drivers_merge_l2');
SW_time = Halekas_SW(:,UT_num);
SW_Bx = Halekas_SW(:,Bx_num);
SW_Bz = Halekas_SW(:,Bz_num);
SW_By = Halekas_SW(:,By_num);
SW_v = Halekas_SW(:,v_num);
SW_n_p = Halekas_SW(:,n_p_num);
SW_dyn_press = 1.6726e-6 *SW_n_p .* SW_v.^2;

end