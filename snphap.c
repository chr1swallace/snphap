#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "cline.h"
#include "snphap.h"


int n_loci = 0;  /* Total number of loci */
int n_phase;    /* The number currently phased */
FILE *warnings; /* File for warning messages */
long n_warn = 0;/* Count of warnings */

int main(int argc, char** argv) {


  /* Additional parameters */

  double memory = 0.0; /* Maximum memory to be consumed by task, in Mbyte*/
  double min_posterior = 0.001; /* Minimum posterior probability */
  double min_prior = 0.0;  /* Minimum prior probability */
  int kill = 5;            /* Cull every ... loci */
  double tol = 0.001;      /* Log likelihood tolerance */
  int maxit = 100;         /* Maximum EM iterations */
  int num = 0;             /* Force numeric output */
  int ss = 0;              /* Tab delimited ouput for spreadsheet */ 
  int wide = 0;            /* Wide output on assignment file */
  int rs = 1;              /* Random starting point */
  long seed = 0;           /* Random number seed */
  int rano = 0;            /* Add loci in random order */
  int revo = 0;            /* Add loci in reverse order */
  int cmpo = 0;            /* Add loci in order of completeness */
  int mafo = 0;            /* Add loci in decreasing MAF order */
  int mmle = 0;            /* Multiple maximizations at last step */
  int mimp = 0;            /* Multiple random imputations */
  int mcmc = 50;           /* Number of MCMC steps between samples */
  int nmh = 0;             /* Locus name header line */
  double dfs = 0.1;        /* Starting df parameter for Dirichlet (*N) */
  double dfe = 0.0;        /* Finishing df parameter for Dirichlet (*N) */
  int ranp = 0;            /* Restart from random prior probabilities */
  int quiet = 0;           /* Suppress screen progress report */
  double of_max = 1.0;     /* Posterior threshold for listing */

  /* Not currently setable */

  int maxs = 10000;         /* Maximum subjects (doesn't matter yet) */
  double head_room = 0.5;  /* Dynamic memory headroom in Mb */

  /* Work variables */

  int n_subject; /* Number of subjects */
  int mpl, mps, mph; /* Memory per locus, per subject, and per haplotype */
  long n_hap, hap_n; /* Number of haplotype records */
  long max_haps; /* maximum number of haplotypes */
  HAP **so_list; /* Subject order haplotype list */
  HAP **ho_list; /* Haplotype order haplotype list */
  HAP **unique;  /* List of unique haplotypes */
  HAP *h1, *h2;
  HAP **h;
  int *order = (int *) 0;    /* Order of adding loci */
  FILE *infile, *outfile;
  char ifname[MAX_FILENAME_LEN], of1name[MAX_FILENAME_LEN],
    of2name[MAX_FILENAME_LEN], tempname[MAX_FILENAME_LEN],
    lfname[MAX_FILENAME_LEN];
  char lname[MAX_NAME_LEN], sfmt[8], c, *cl, *cname;
  char **names;
  CODE *coding;
  int  res, j, k, iter, carry_on, nlen, mi1, mi2, mi3, ladd;  
  long i;
  time_t now;
  double memmax, logl, lastl, best_logl, df, df_step; 
  double *p_unique;

  /* Get parameters, calculate limits, open files */

  *lfname = (char) 0;
  if (
      get_flag(argc, argv, "l", 1, &n_loci) < 0 ||
      get_flag(argc, argv, "mb",3 , &memory) < 0 ||
      get_flag(argc, argv, "k", 1, &kill) < 0 ||
      get_flag(argc, argv, "pr",3 , &min_prior) < 0 ||
      get_flag(argc, argv, "po",3 , &min_posterior) < 0 ||
      get_flag(argc, argv, "to",3 , &tol) < 0 ||
      get_flag(argc, argv, "th",3 , &of_max) < 0 ||
      get_flag(argc, argv, "i",1 , &maxit) < 0 ||
      get_flag(argc, argv, "nu",0 , &num) < 0 ||
      get_flag(argc, argv, "nh",0 , &nmh) < 0 ||
      get_flag(argc, argv, "nf",4, lfname) < 0 ||
      get_flag(argc, argv, "ss",0 , &ss) < 0 ||
      get_flag(argc, argv, "w", 0, &wide) < 0 ||
      get_flag(argc, argv, "rs",0 , &rs) < 0 ||
      get_flag(argc, argv, "rp",0 , &ranp) < 0 ||
      get_flag(argc, argv, "ro",0 , &rano) < 0 ||
      get_flag(argc, argv, "rv",0 , &revo) < 0 ||
      get_flag(argc, argv, "co",0 , &cmpo) < 0 ||
      get_flag(argc, argv, "mo",0 , &mafo) < 0 ||
      get_flag(argc, argv, "sd",5 , &seed) < 0 ||
      get_flag(argc, argv, "mm",1 , &mmle) < 0 ||
      get_flag(argc, argv, "mi",1 , &mimp) < 0 ||
      get_flag(argc, argv, "q", 0, &quiet) < 0 ||
      (mi1 = get_flag(argc, argv, "mc", 1, &mcmc)) <0 ||
      (mi2 = get_flag(argc, argv, "ds", 3, &dfs)) <0 ||
      (mi3 = get_flag(argc, argv, "de", 3, &dfe)) <0  ||
      !get_arg(argc, argv, ifname)
      ) {
    fprintf(stderr, "\nSNPHAP, Version ");
    fprintf(stderr, VERSION);
    fprintf(stderr, "\nUse:\n\n");
    fprintf(stderr, "%s", argv[0]);
    fprintf(stderr, "\t[-i # -k # -l # -mb # -mi # -mm # -nu -nh -pr #  -po # \n");
    fprintf(stderr, "\t -q -[co|mo|ro|rv] -rs -sd # -ss -w -to # -th # -mc # -ds # -de #]");
    fprintf(stderr, "\n\t\t\tinput-file [output-file output-file]\n");
    fprintf(stderr, "\nWhere flag settings (and defaults) are:\n\n");
    fprintf(stderr, "\t-i\tMaximum EM iterations (%d)\n", 
	    maxit);
    fprintf(stderr, "\t-k\tCull rare haplotype assignments every k loci (%d)\n",kill);
    fprintf(stderr, "\t-l\tNumber of loci (if no header line in file)\n");
    fprintf(stderr, "\t-mb\tMaximum dynamic storage to be allocated, in Mb (%.1lf)\n", 
	    memory);
    fprintf(stderr, "\t-mi\tCreate multiple imputed datasets (%d).  ", mimp);
    fprintf(stderr, "If set >0:\n");
    fprintf(stderr, "\t  -mc\tNumber of MCMC steps between samples (%d)\n",
	    mcmc);
    fprintf(stderr, "\t  -ds\tStarting value of Dirichlet prior parameter (%.2lf*N)\n",
	    dfs);
    fprintf(stderr, "\t  -de\tFinishing value of Dirichlet prior parameter (%.2lf*N)\n",
	    dfe);
    fprintf(stderr, "\t-mm\tRepeat final maximization multiple times (%d)\n",
	    mmle);
    fprintf(stderr, "\t-nu\tForce numeric coding of alleles in output\n");
    fprintf(stderr, "\t-nh\tLocus names form first line of input\n");
    fprintf(stderr, "\t-pr\tPrior (ie population) probability threshold (%0.2g)\n", 
	    min_prior);
    fprintf(stderr, "\t-po\tPosterior probability threshold (%0.2g)\n", 
	    min_posterior);
    fprintf(stderr, "\t-q\tQuiet operation (off)\n");
    fprintf(stderr, "\t-co\tLoci added in order of completeness (off)\n");
    fprintf(stderr, "\t-mo\tLoci added in MAF orderr (off)\n");
    fprintf(stderr, "\t-ro\tLoci added in random order (off)\n");
    fprintf(stderr, "\t-rv\tLoci added in reverse order (off)\n");
    fprintf(stderr, "\t-rs\tRandom starting points for each EM iteration (off)\n");
    fprintf(stderr, "\t-sd\tSet seed for random number generator (use date+time)\n");
    fprintf(stderr, "\t-ss\tTab-delimited speadsheet file output (off)\n");
    fprintf(stderr, "\t-to\tLog-likelihood convergence tolerance (%0.2g)\n", 
	    tol);
    fprintf(stderr, "\t-th\tPosterior probability threshold for output (%0.4g)\n", 
	    of_max);
    fprintf(stderr, "\t-w\tWide format for output file 2 (haplotype assignments)\n");
    fprintf(stderr, "\n");
    return 1;
  }

  warnings = fopen("snphap_warnings", "w");
  if (!warnings) {
    fprintf(stderr, "Couldn't open warnings file\n");
    return 1;
  }

  /* Optional arguments */

  get_arg(argc, argv, of1name);
  get_arg(argc, argv, of2name);

  /* If both missing give default name to first */

  if (!of1name[0] && !of2name[0]) {
    strcpy(of1name, "snphap.out");
  }
    
  /* Check for unrecognized arguments */

  cl = unrec(argc, argv);
  if (cl) {
    fprintf(stderr, "Unexpected flag or argument: %s\n", cl);
    return 1;
  }

  /* Any incompatible arguments */

  if ((rano+revo+cmpo+mafo)>1) {
    fprintf(stderr, "No more than one of -rv -ro -co -mo flags can be set\n");
    return 1;
  } 
  if (mimp<=0 && (mi1 || mi2 || mi3)) {
    fprintf(stderr, "-mc, -ds, -de parameters are only legal if -mi is set\n");
    return 1;
  }
  if (ranp && mmle<=0) {
    fprintf(stderr, "-rp option only relevant with -mm # option\n");
    return 1;
  }
  
  time(&now);
  printf("\nSNPHAP, Version ");
  printf(VERSION);
  if (now != -1) {
    printf("\nRun on  %s", asctime(localtime(& now)));
    fprintf(warnings, "snphap warnings: %s", asctime(localtime(& now)));
  }
  else {
    printf(": time not available\n");
  };
  memmax = memavail(1024)/1.0e6;
  if (memory<head_room || memory> memmax) {
    printf("Dynamic memory available %.1lf Mbyte\n", memmax);
    memory = memmax;
  }
 
  if (nmh) { 
    if (*lfname) {
      fprintf(stderr, "-nh and -nf cannot be used together\n");
      return 1;
    }
    infile = fopen(ifname, "r");
    if (!infile) {
      fprintf(stderr, "Couldn't open input file\n");
      return 1;
    }
    printf("\nLocus names read from input file:");
    n_loci = field_count(infile); 
    names = calloc(n_loci, sizeof(char *));
    if (!names) 
      goto no_room;
    nlen = 0;
    for (j=0; j<n_loci; j++) {
      fscanf(infile, " %s", lname);
      cname = malloc(1+strlen(lname));
      if (!cname) goto no_room;
      strcpy(cname, lname);
      if (strlen(lname)>nlen)
	nlen = strlen(lname);
      names[j] = cname;
      printf(" %s", cname);
    }
    printf("\n");
  }   
  else if (*lfname) {
    printf("\nLocus names read from file %s:", lfname);
    infile = fopen(lfname, "r");
    if (!infile) {
      fprintf(stderr, "Couldn't open locus name file\n");
      return 1;
    }
    n_loci = 0;
    while ((c=fgetc(infile))!=EOF) {
      if (c=='\n') 
	n_loci++;
    }
    rewind(infile);
    names = calloc(n_loci, sizeof(char *));
    if (!names) 
      goto no_room;
    sprintf(sfmt, " %%%ds", MAX_NAME_LEN-1);
    nlen = 0;
    for (j = 0; j<n_loci; j++) {
      fscanf(infile, sfmt, lname);
      cname = malloc(1+strlen(lname));
      if (!cname) goto no_room;
      strcpy(cname, lname);
      names[j] = cname;
      printf(" %s", cname);
      if (strlen(lname)>nlen)
	nlen = strlen(lname);
      while(fgetc(infile)!='\n');
    }
    printf("\n");
    fclose(infile);
    infile = fopen(ifname, "r");
    if (!infile) {
      fprintf(stderr, "Couldn't open input file\n");
      return 1;
    }
  }
  else {
    infile = fopen(ifname, "r");
    if (!infile) {
      fprintf(stderr, "Couldn't open input file\n");
      return 1;
    }
    if (!n_loci) {
      n_loci = field_count(infile) - 1;
      if (n_loci % 2) {
	fprintf(stderr, "Error in format of input file\n");
	fprintf(stderr, "(Should be an id plus two fields for each locus)\n");
	return 1;
      }
      n_loci /= 2;
    }
    names = (char **) 0;
    nlen = 1;
  }
  printf("\nNumber of loci                  : %10d\n", n_loci);
  mpl = sizeof(CODE);
  if (names)
    mpl += sizeof(char *) + nlen;
  if (revo || rano || cmpo) 
    mpl += 2*sizeof(int);
  if (mafo)
    mpl += (sizeof(int)+2*sizeof(double));

  /* Dynamic memory calculations */

  mps = 0;
  mph = 3*sizeof(HAP *) + sizeof(HAP); 
  if (mmle > 0 || mimp > 0) 
    mph += sizeof(double);
  max_haps = ((memory-head_room)*1.e6 - n_loci*mpl - maxs*mps)/mph;
  printf("Maximum memory usage            : %10.1lf Mb\n", memory);
  printf("Maximum haplotype instances     : %10ld\n", max_haps);
  if (rano || rs || mimp>0) {
    if (!seed) {
      seed = (long) now;
    }
    printf("Seed for random numbers         : %10ld\n", seed);
    SEED (seed);
  }
  if (rs) {
    rs = 1;
    printf("\nEach EM iteration will be started from a random imputation\n");
  }

  /* set up main arrays */

  so_list = calloc(max_haps, sizeof(HAP*));
  if (!so_list) 
    goto no_room;
  ho_list = calloc(max_haps, sizeof(HAP*));
  if (!ho_list) 
    goto no_room;
  coding = (CODE *) calloc(n_loci, sizeof(CODE));
  if (!coding) 
    goto no_room;
  else {
    for (j=0; j<n_loci; j++) {
      coding[j].anum = 0;
      coding[j].one = coding[j].two = (char) 0;
    }
  }
 
  /* Read in data */

  n_subject = n_hap = 0;
  h = so_list;
  while ( (res = gt_read(infile, coding, &h1, &h2)) ) {
    switch (res) {
    case 1:
      n_subject++;
      *h++ = h1;
      *h++ = h2;
      n_hap += 2;
      break;
    case 2:
      fprintf(stderr, "Skipping this subject(%s)\n", h1->id);
      break;
    case 3:
      fprintf(stderr, "Exiting on data input error\n");
      return 1;
    case 4:
      printf("End of file assumed after %d subjects\n", n_subject);
    }
    if (res==4) break;
  }
  printf("Subjects read from data file    : %10d\n", n_subject);

  /* Reorder if necessary */

  if (rano) {
    order = (int *) calloc(n_loci, sizeof(int));
    if (order) {
      ranord(n_loci, order);
      /*
      printf("Loci added in the order         :");
      for (i=0; i<n_loci; i++)
	printf(" %d", 1+order[i]);
      printf("\n");
      */
      printf("Loci added in random order");
    }
    else {
      goto no_room;
    }
    reorder(n_hap, so_list, order);
  }
  else if (revo) {
    order = (int *) calloc(n_loci, sizeof(int));
    printf("Loci will be added in reverse oder\n");
    if (order) {
      for (i=0; i<n_loci; i++) 
      order[i] = n_loci-i-1;
    }
    else {
      goto no_room;
    }
    reorder(n_hap, so_list, order);
  }
  else if (cmpo) {
    order = (int *) calloc(n_loci, sizeof(int));
    printf("Loci will be added in order of data completeness\n");
    if (order) {
      compord(n_hap, so_list, order);
    }
    else {
      goto no_room;
    }
    reorder(n_hap, so_list, order);
  }
  else if (mafo) {
    order = (int *) calloc(n_loci, sizeof(int));
    printf("Loci will be added in decreasing order of MAF\n");
    if (order) {
      maford(n_hap, so_list, order);
    }
    else {
      goto no_room;
    }
    reorder(n_hap, so_list, order);
  }

  /* Show thresholds and convergence criteria */

  if (kill<=n_loci) {
    printf("\nRare haplotype assignments are culled every %d loci\n", kill);
    printf("Thresholds for trimming\n");
    printf("\tPrior probability       : %10.2lg\n", min_prior);
    printf("\tPosterior probability   : %10.2lg\n", min_posterior);
  }
  printf("\nEM convergence criteria\n");
  printf("\tMaximum iterations      : %10d\n", maxit);
  printf("\tTolerated LLH change    : %10.2lg\n", tol);
  if (mmle > 0) {
    printf("Repetions of final EM iteration : %10d\n", mmle);
  }
  if (mimp > 0) {
    printf("Multiple imputed datasets       : %10d\n", mimp);
    printf("\tMCMC steps              : %10d\n", mcmc);
    dfs *= (2*n_subject);
    dfe *= (2*n_subject);
    df_step = (dfs - dfe)/(double) mcmc;
    printf("\tStarting Dirichlet df   : %10.2lf\n", dfs);
    printf("\tFinishing Dirichlet df  : %10.2lf\n", dfe);
  }
  if (of1name[0]) {
    printf("Haplotype frequencies output to : %10s", of1name);
    if (mimp>0) 
      printf("(.*)\n");
    else
      printf("\n");
  }
  if (of2name[0]) {
    printf("Haplotype assignments output to : %10s", of2name);
    if (mimp>0) 
      printf("(.*)\n");
    else
      printf("\n");
  }

  /* Introduce loci to resolve one at a time */

  if (!quiet) {
    printf("\nProgress of maximum likelihood estimation:\n");
    printf("\nLocus added\tLoci Phased\tht instances\tIteration\tLog Likelihood");
    printf("\n-----------\t-----------\t------------\t---------\t--------------\n");
  }
  n_phase = 0;
  while (n_phase < n_loci) {
    ladd = order? order[n_phase]+1: n_phase+1;
    
    /* Expand by  possible phase at last locus */

    n_hap = hap_expand(n_hap, max_haps, so_list, rs);
    if (!n_hap) goto no_room;
    if (!quiet) {
      printf("\r%11d\t%11d\t%12ld\t%9s", ladd, n_phase, n_hap, "(sorting)");
      fflush(stdout);
    }

    /* Create list sorted by haplotype */
    
    for (i=0; i<n_hap; i++) 
      ho_list[i] = so_list[i];
    qsort(ho_list, n_hap, sizeof(HAP *), cmp_hap);

    /* Iteration to fit haplotypes as far as locus n_phase */

    carry_on =  1;
    iter = 0;
    while (carry_on) {
      iter++;
      hap_prior(n_hap, ho_list, min_prior);
      hap_n = hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
      carry_on = (iter==1) || ((logl - lastl)>tol);
      if (carry_on && (iter==maxit)) {
	carry_on = 0;
	fprintf(warnings, "No convergence in %d iterations", iter);
	fprintf(warnings, "\t(at %d-locus step)\n", 
		n_phase);
	n_warn++;
      }
      lastl = logl;
      if (!quiet) {
	printf("\r%11d\t%11d\t%12ld\tEM%7d\t%14.3lf", ladd, n_phase, n_hap, 
	       iter, logl);
	fflush(stdout);
      }
    }
    n_hap = hap_posterior(n_hap, so_list, min_posterior, &logl,
			  !((n_phase+1) % kill) );
  }
  if (!quiet) { 
    printf("\r%11d\t%11d\t%12ld", ladd, n_phase, n_hap);
    fflush(stdout);
  }

  /* Wrap up after last iteration and compute unique haplotypes */

  for (i=0; i<n_hap; i++) 
    ho_list[i] = so_list[i];
  qsort(ho_list, n_hap, sizeof(HAP *), cmp_hap);
  hap_prior(n_hap, ho_list, min_prior);
  hap_n = n_unique_haps(n_hap, ho_list);
  unique = (HAP **) calloc(hap_n, sizeof(HAP *));
  if (!unique)
    goto no_room;
  unique_haps(n_hap, ho_list, unique);
  if (mmle>0 || mimp>0) {
    p_unique = (double *) calloc(hap_n, sizeof(double));
    if (!p_unique)
      goto no_room;
    for (i=0; i<hap_n; i++) {
      p_unique[i] = unique[i]->prior;
    }
  }

  /* Repeat last iteration, saving most likely solution */

  if (mmle > 0) {
    best_logl = logl;
    for (j=0; j<mmle; j++) {
      carry_on =  1;
      iter = 0;
      if (ranp) {
	hap_prior_restart(n_hap, ho_list);
	hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
      }
      else {
	hap_posterior_restart(n_hap, so_list);
      }
      while (carry_on) {
	iter++;
	hap_prior(n_hap, ho_list, min_prior);
	hap_posterior(n_hap, so_list, min_posterior, &logl, 0);
	carry_on = (iter==1) || ((logl - lastl)>tol);
	if (carry_on && (iter==maxit)) {
	  carry_on = 0;
	  fprintf(warnings, "No convergence in %d iterations", iter);
	  fprintf(warnings, "\t(at %d/%04d-locus step)\n", 
		  n_loci, j+1);
	  n_warn++;
	}
	lastl = logl;
	if (!quiet) { 
	  printf("\r%11d\t%6d/%04d\t%12ld\tEM%7d\t%14.3lf", ladd, n_loci, j+1, 
	       n_hap, iter, logl);
	  fflush(stdout);
	}
      }
      if (logl > best_logl) {
	best_logl = logl;
	for (i=0; i<hap_n; i++) {
	  p_unique[i] = unique[i]->prior;
	}
	if (!quiet) {
	  printf("\r\t\t\t\t\t\t\t\t\t(%.3lf)", best_logl);
	  fflush(stdout);
	}
      }	
    }
    hap_prior_restore(n_hap, ho_list, p_unique);
    if (quiet)
      printf("\n\nBest solution has log likelihood %.3lf", best_logl);
  }

  qsort(unique, hap_n, sizeof(HAP *), more_probable);

  /* Output MLE etc */

  if (!quiet) 
    printf("\n");
  printf("\n");
  if (of1name[0]) {
    outfile = fopen(of1name, "w");
    if (!outfile) 
      goto open_error;
    if (!ss) {
      fprintf(outfile, "snphap listing: %s", asctime(localtime(& now)));
    }
    j = hap_write(outfile, n_loci, names, coding, order, hap_n, unique, 
		  0, 0.0, num, ss);
    fclose(outfile);
    printf("%d haplotypes written to output file %s\n", j, of1name);
  }
  if (of2name[0]) {
    outfile = fopen(of2name, "w");
    if (!outfile) 
      goto open_error;
    j = hap_write(outfile, n_loci, names, coding, order, n_hap, so_list, 
		  wide+1, of_max, num, ss);
    fclose(outfile);
    printf("%d possible assignments to %d subjects written to output file %s\n",
	   j, n_subject, of2name);
  }

  /* Multiple imputation by IP algorithm */

  if (mimp > 0 && !quiet) {
    printf("\nProgress of multiple imputation:\n");
    printf("\nImputation\tIP step\t    Prior df");
    printf("\n----------\t-------\t    --------\n");
  }
  for (j=1; j<=mimp; j++) {
    hap_prior_restore(n_hap, ho_list, p_unique);
    df = dfs;
    for (iter=0; iter<mcmc; iter++) {
      df -= df_step;
      if (!quiet) {
	printf("\r%11d\t%7d\t%12.2lf", j, iter+1, df);
	fflush(stdout);
      } 
      sample_posterior(n_hap, so_list); /* I step */
      sample_prior(n_hap, ho_list, df); /* P step */
    }
    if (of1name[0]) {
      sprintf(tempname,"%s.%03d", of1name, j);
      outfile = fopen(tempname, "w");
      hap_write(outfile, n_loci, names, coding, order, hap_n, unique, 
		  0, 0.0, num, ss);
      fclose(outfile);
    }
    if (of2name[0]) {
      sprintf(tempname,"%s.%03d", of2name, j);
      outfile = fopen(tempname, "w");
      hap_write(outfile, n_loci, names, coding, order, n_hap, so_list, 
		  wide+1, 0.0, num, ss);
      fclose(outfile);
    }    
  }
  if (!quiet) 
    printf("\n");

  if (mimp>0 && of1name[0]) {
    printf("\nSamples from posterior distribution of haplotype frequencies ");
    printf("written to \n\tfiles %s.001 ... %s.%03d\n", of1name, of1name, 
	   mimp);
  }
  if (mimp>0 && of2name[0]) {
    printf("Multiply imputed datasets written to files %s.001 ... %s.%03d\n",
	   of2name, of2name, mimp);
  }

  /* Describe numerical recoding */

  if (num) {
    printf("\nNumerical recoding of alleles was forced:\n");
    for (j=k=0; j<n_loci; j++) {
      if (coding[j].anum == 2) {
	k++;
	printf("%s: \t", names[j]);
	c = allele_code(1, coding[j]);
	if (c)
	  printf("1=%c", c);
	c = allele_code(2, coding[j]);
	if (c)
	  printf(", 2=%c", c);
	printf("\n");
      }
    }
    if (!k) 
      printf("\t... but no recoding was necessary\n");
  }

  /* Normal return */

  if (n_warn) 
    printf("\nNote: %ld messages were written to the warnings file\n", n_warn);
  else 
    printf("\n");
  return 0;

  /* Error conditions */

 no_room:
  fprintf(stderr, "Insufficient memory\n");
  return 1;

 open_error:
  fprintf(stderr, "Error opening output file\n");
  return 1;

}

  

  
