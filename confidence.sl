% List of paramers to run on.
variable p = [1,2,3,4,5,6,7,8,11,15,16,17,19];

variable all_good, good, i, lo, hi;

do {
	vmessage("\n\n--Starting to determine confidence intervals--\n");
	all_good = 1;
	fp = fopen("Conf_2009.csv","w");% Create the file
	for (i=0; i<length(p); i++) {
		do {
			good = 1;
			(lo, hi) = conf(p[i]);
			vmessage("  ... parameter %i = %e (%e - %e)", p[i], get_par(p[i]), lo, hi);
			() = fprintf(fp, "%s,%e,+%e,-%e\n", get_par_info(p[i]).name, get_par(p[i]), hi-get_par(p[i]), get_par(p[i])-lo);
			if (lo == hi) {
				% this means that a better fit was found
				vmessage("Found a better fit!");
				() = fit_counts;
				good = 0;
				all_good = 0;
			}
		} while (good == 0);
		vmessage("Parameter %i = %e (%e - %e)", p[i], get_par(p[i]), lo, hi);
	}
	() = fclose(fp);
} while (all_good == 0);

save_par("conf_2009.par");
()=eval_counts;
list_par;
