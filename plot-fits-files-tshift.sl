% Load required data
tm1 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0000/lag.fits");
tm2 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0020/lag.fits");
tm3 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0040/lag.fits");
tm4 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0050/lag.fits");
tm5 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0080/lag.fits");
tm6 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0100/lag.fits");
tm7 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0125/lag.fits");
tm8 = fits_read_table ("lag-files/i53/a0.998/h02.0+03.0/g2.7+2.3/A01.0/x2.0/p00.0/b02.0/q0.5/tmax1000/tshift0150/lag.fits");


% Plot results
%open_plot("lag-files-delta-tshift2.png/png");
%variable pid = open_plot("lag-files-delta-tshift2.ps/cps");
charsize(1.2); %good font size
nice_width; %good line ratios
_pgscf(2); %make font roman
xlog;
ylin;
popt.res = 1;
charsize(1.2);
ylabel("Time Lag (t\dg)");
xlabel("Frequency (1/t\dg)");
%title("Frequency dependent Lags (t_max = 1000, t_shift = 0 - 40)");
connect_points(1);
plotxy(tm1.fq,tm1.lag;dcol=2,dsym=1, xrange=[1e-5,1.0], yrange=[-50, 50]);
plotxy(tm2.fq,tm2.lag;dcol=3,dsym=1,oplt=1);
plotxy(tm3.fq,tm3.lag;dcol=4,dsym=1,oplt=1);
plotxy(tm4.fq,tm4.lag;dcol=5,dsym=1,oplt=1);
plotxy(tm5.fq,tm5.lag;dcol=6,dsym=1,oplt=1);
plotxy(tm6.fq,tm6.lag;dcol=7,dsym=1,oplt=1);
plotxy(tm7.fq,tm7.lag;dcol=8,dsym=1,oplt=1);
plotxy(tm8.fq,tm8.lag;dcol=1,dsym=1,oplt=1);
nice_width;
_pgmove(-10000,0);
_pgsls(2);
_pgslw(1);
_pgdraw(1,0);
%close_plot;
