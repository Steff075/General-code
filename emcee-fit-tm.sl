% Two blobs model fitting routines for EMCEE following the Mike Nowak REF: https://www.sternwarte.uni-erlangen.de/wiki/index.php/Emcee
% This script includes the fitting routine (written by Andy Young) - fitting the 2 blobs model to a given data set via the best parameter fits found on BlueCrystal
%
% Date: December 2020
%
% Should we run the fit? Non-zero => run fit
variable run_fit = 1;

% Filename prefix for observation identifier
variable obsid = "1H0707_comb_fmax";

% Some global variables and paths (update to your own path)
% variable path_to_lags = "/mnt/storage/scratch/phajy/x-ray-reverberation";
% Assumes the lags are in a subdirectory of the current directory
%variable path_to_lags = "/mnt/storage/scratch/so14935/x-ray-rev4/x-ray-reverberation";
variable path_to_lags = "/data/typhon1/reverb";
variable path_to_fits = getcwd;

% load emcee asuming that "isis_emcee.sl" is located in the cwd
require("isis_emcee");
require("isisscripts");

% Debugging - non-zero to output more information
% higher value means more diagnostic information
variable debug = 0;

% Only need to cache the last file read in because the script steps through
% discrete parameter values. The number of parameter combinations has between
% reduced to make this problem tractable. This vastly simplifies the code from
% previous versions and reduces the memory requirements significantly.

variable last_file_used = "";
variable last_lag_used = NULL;

% Lag-frequency fitting function
define lag_freq_fit(lo, hi, par)
{
	variable norm = par[0];
	variable mass = 10.0^par[1];

	variable i = par[2];
	variable a = par[3];
	variable h1 = par[4];
	variable h2 = par[5];
	variable g1 = par[6];
	variable g2 = par[7];
	variable abund = par[8];
	variable xi = par[9];
	variable p = par[10];
	variable b = par[11];
	variable q = par[12];
	variable tmax = par[13];
	variable tshift = par[14];

	% Now we need to generate the appropriate filename and see if it exists
	variable fp;
	variable filename = sprintf("lag-files/i%02.0f/a%05.3f/h%04.1f+%04.1f/g%03.1f+%03.1f/A%04.1f/x%03.1f/p%04.1f/b%04.1f/q%03.1f/tmax%04.0f/tshift%04.0f/lag.fits", i, a, h1, h2, g1, g2, abund, xi, p, b, q, tmax, tshift);
	if (debug>5) vmessage("Using filename %s", filename);

	variable t_g = 492.703806 * (mass / 1.0E8);

	% Read in the lag as a function of frequency
	% Note that units are tg^-1 (freq) and tg (lag)

	variable spec;
	variable ft;
	variable zero = 0;

	% see if we've already loaded in this file
	if (last_file_used == filename) {
		% we have already seen this file
		% retrieve the saved values from memory
		if (debug>5) vmessage("Using cached lag");
		ft = @last_lag_used;
	} else {
		% we need to read in the file
		% check if file exists
		if (debug>1) vmessage("checking " + path_to_lags + "/" + filename);
		fp = fopen(path_to_lags + "/" + filename, "r");
	  if (fp == NULL)
	  {
	  	% file does not exist
	    if (debug>0) vmessage("File %s does not exist so returning zero", filename);
			last_file_used = "";
			last_lag_used = NULL;
	   } else {
	   	() = fclose(fp);

    	if (debug > 0) vmessage("Reading file %s", filename);
			ft = fits_read_table(path_to_lags + "/" + filename);
			last_file_used = filename;
			last_lag_used = @ft;

			% Check for NaN which indicates a problem
			variable cn = where(isnan(ft.lag) == 1);
			if (length(cn)>0) {
				vmessage("File %s has NaN in it (%i).", filename, length(cn));
				plot_pause;
			} else {
				if (debug > 5) vmessage("Files %s is good.", filename);
			}
		}
	}

	if (last_file_used != "") {
		% Shift all of the frequency bins
		variable fq = ft.fq / t_g;

		% Scale all time lags
		variable lag = ft.lag * t_g;

		% Overall normalizaton (which accounts for dilution, etc.)
		lag = lag * norm;

		% Frequency should now be in units of Hz
		% Lag should now be in units of seconds

		spec = rebin(lo, hi, _A(fq), _A(make_hi_grid(fq)), reverse(lag));

		% Check for NaN which indicates a problem
		cn = where(isnan(spec)==1);
		if (length(cn)>0) {
			vmessage("Spectrum from file %s has NaN in it.", filename);
			plot_pause;
		}

	} else {
		spec = Double_Type[length(lo)];
	}

	return(spec);
}

% New version will be:
add_slang_function("lag_freq", ["norm", "log_mass", "i", "a", "h1", "h2", "Gamma1", "Gamma2", "Fe", "xi", "p", "b", "q2", "t_max", "t_shift"]);

% Load the data.
Ignore_PHA_Response_Keywords = 1;
variable d = load_data("/data/typhon1/reverb/1H0707-495/1H0707-495-time-lag-fits/1H0707-495-comb-time-lags-fmax.fits");
variable r = load_rmf("/data/typhon1/reverb/1H0707-495/1H0707-495-time-lag-fits/1H0707-495-comb-time-lags-fmax.rsp");
assign_rmf(r, d);

%load_par("1H0707_0511580201_A_best_fit.par");
%load_par("conf_comb.par");
load_par("conf_comb_test.par");
Label_By_Default = 0;
popt.xlabel="Frequency (Hz)";
popt.ylabel="Time lag (s)";
popt.yrange=[-50, 200];
popt.xrange=[6e-5, 0.01];
%xlog;

% Instead of using fit_fun, run the emcee sampler for 500 walkers and 1000 iterations
%emcee(200, 50000; output="1H0707_comb_TM_fit_chain.fits", serial);

% Read the chain-file
(nw, nfreepar, freepar, chain_statistic, chain) = read_chain("1H0707_comb_TM_fit_chain.fits");

close_plot; % start with a fresh plot window
_for i (0, nw-1, 1) { % loop over all walkers and plot their path
  chain_plot(1,i, nw, chain; oplot); 
};

%Try to get better labelsgedit 
%chain_plqot(1,i, nw, chain; oplot);
%charsize(1.4);
%nice_width;


% Make good image of the paths of the walkers
variable pid = open_plot("1H0707_comb_TM_fit_walkers.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
nice_width;
chain_plot(1,i,nw,chain);
close_plot(pid);


% Obtain acceptance Rate 
frac_update = read_chain("1H0707_comb_TM_fit_chain.fits"; frac_update);
xlabel("Iteration Step");
ylabel("Acceptance Rate");
popt.dsym = [6];
popt.dcol = 4; 
plot([0:length(frac_update)-1], frac_update);
plot(([0:length(frac_update)-1], frac_update); dsym=6, dcol=2);

% Make good image of the acceptance Rate
variable pid = open_plot("1H0707_comb_TM_fit_mass_accept_rate.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
xlabel("Iteration Step");
ylabel("Acceptance Rate");
plot([0:length(frac_update)-1], frac_update);
close_plot(pid);

% Obtain the probability distributions of the norm
chain_hist(0, chain; fcut = .3, pmin = 0, pmax = 0.0001, nbin = 200); % norm

% Make good image of the Norm probability
variable pid = open_plot("1H0707_comb_TM_fit_norm_prob.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
chain_hist(0,chain; fcut = 0.3,min = 0, pmax = 0.0001, nbin=200);
close_plot(pid);

% Obtain the probability distributions of the log_mass
chain_hist(1, chain; fcut = 0.3, pmin = 6.0, pmax = 6.5, nbin = 400); % log_mass

% Make good image of the Mass probability
variable pid = open_plot("1H0707_comb_TM_fit_mass_prob.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
chain_hist(1,chain; fcut = .3, pmin = 6.0, pmax = 6.5, nbin = 200);
close_plot(pid);

% Contours: paths in 2D 
%where 1 = x and 0 = y axis parameters
chain_hist2d(1, 0, chain; nbinx = 200, nbiny = 200, xlo = 6.0, xhi = 6.5, ylo = 1e-5, yhi = 4e-5);

%Make good image of the Contours
variable pid = open_plot("1H0707_comb_TM_fit_chain_hist2d.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
chain_hist2d(1,0,chain; nbinx = 200, nbiny = 200, xlo = 6.0, xhi = 6.5, ylo = 1e-5, yhi = 4e-5);
close_plot(pid);

(frac_update, min_statistic, med_statistic, max_statistic) = read_chain("1H0707_comb_fmax_chain.fits"; chain_stats);

% Chain statistics starting at s=iteration where the fit converged
s=1000; % start at this iteration step
min(min_statistic[[s:]]);

% Plot the statistic array
n=length(max_statistic)-1;
yrange(min(min_statistic[[s:]])*.9, max(max_statistic[[s:]])*1.1);
variable pid = open_plot("1H0707_comb_TM_fit_statistics.ps/cps");
charsize(1.3); %good font size
_pgscf(2); %make font roman
plot([s:n], min_statistic[[1000:]]);	
oplot([s:n], max_statistic[[1000:]]);
oplot([s:n], med_statistic[[1000:]]);
%close_plot(pid);


% Error analysis

variable p = 1; % index of the parameter of interest (here: norm(0))
variable lvl = 0.90; % confidence interval to find
variable fcut = 0.2; % throw away the first fcut-fraction of the chain
variable nbins = 400; % number of histogram bins
chain_hist(p, chain; fcut = fcut, nbin = nbins);
%close_plot;

variable pmin, pmax;
(pmin,pmax,,) = cursor_box();

pmin;
pmax;

% histogram to get the parameter's probability distribution

variable lo, hi;
(lo, hi) = linear_grid(pmin, pmax, nbins);
variable hist = 1.*histogram(chain[p], lo, hi);

% find maximum
variable nmax = where_max(hist)[0];

% do a linear regression of the first and last n bins
variable n = 10;
variable ind = [[0:n-1], [nbins-n:nbins-1]];

% linear regression
variable lin = linear_regression(0.5*(lo+hi)[ind], hist[ind]);

% subtract the linear continuum from the probability distribution
variable newhist = hist - (lin.a + lin.b*.5*(lo+hi));

% finally, normalize!
newhist /= sum(newhist);

% Determine confidence interval
variable li = [:nmax], ri = [nmax+1:];
variable conf_min = lo[li][wherefirst(cumsum(newhist[li]) > .5 - .5*lvl)];
variable conf_max = hi[ri][wherefirst(cumsum(newhist[ri]) > .5*lvl)];
variable conf_mean = conf_min + ((conf_max - conf_min) / 2);
conf_min;
conf_max;
conf_mean;

% Plot the parameter with conf regions

%variable pid = open_plot("1H0707_comb_fmax_logmass_conf68.ps/cps");
variable pid = open_plot("1H0707_comb_TM_fit_logmass_conf90.ps/cps");
chain_hist(1, chain; fcut = fcut, nbin = nbins, pmin=6.00, pmax=6.50, xlabel='logM');
_pgscf(2); %make font roman
charsize(1.4); %good font size
nice_width;
xlabel="logM";
ylabel="Probability";
_pgmove(conf_min, 0);
_pgsls(2); % line style
_pgsci(1); % line color
_pgslw(1); % line width
_pgdraw(conf_min, 100);

_pgmove(conf_max, 0);
_pgsls(2);
_pgsci(1);
_pgslw(1);
_pgdraw(conf_max, 100);

_pgmove(conf_mean, 0);
_pgsls(1);
_pgsci(4);
_pgslw(2);
_pgdraw(conf_mean, 100);

%close_plot;
close_plot(pid);



 
