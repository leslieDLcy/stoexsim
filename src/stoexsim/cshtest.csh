#!/bin/tcsh

# echo $pwd
sed -i -e s/"#SITEAMP#"/52/g make_ds_aleatoric/crustal_amps.txt
# echo 'test begins'






# ##### basic math operation #####

# # # a test of log and exponential operations
# # set a_num=10

# # # this a test log10 operation
# # set log_test=`echo "" | awk '{print ((log('$a_num')/log(10.)))}'`
# # # set log_sd_mean=`echo "" | awk '{print (('$sd'))}'`

# # echo $log_test

# # # now let's try the 'exp' operation
# # set exp_a_num=`echo "" | awk '{print ((10^('$log_test')))}'`
# # echo $exp_a_num






# ##### Draw from a uniform distribution #####

# # # fprintf(stderr, "\n\nusage: ndev <mean> <stdev> $RANDOM [min] [max]\n\n");
# # set RANDOM=`bash -c 'echo $RANDOM'`	# select a psuedo-random number
# # set num=`./random/unidev 1 1 $RANDOM -10 10` #normal distibution

# # python -c "import numpy; print(numpy.random.randint(5,10))"

# set num=`python -c "import numpy; print(numpy.random.uniform(0.1, 0.5))"`
# echo "from python we generate:"$num




# # ##### Draw a sample of `SD` from a PDF
# # set DYN_RAN=`bash -c 'echo $RANDOM'` 
# # # echo $DYN_RAN

# # set RANDOM=`bash -c 'echo $RANDOM'`	# select a psuedo-random number

# # # sd value in natural scale
# # set sd=50 

# # # change sd into log10 scale as the mean
# # set sd_log10=`echo "" | awk '{print (( log('$sd')/log(10.) ))}'`
# # echo "log10 value of 'sd=50' ="$sd_log10

# # # mannually set up the std for the lognormal distribution
# # set sd_dev=0.31

# # # Hint: -Normal distribution - usage: ndev <mean> <stdev> $RANDOM [niter] [min max]\n
# # set sd_sample_log10=`./random/normdev $sd_log10 $sd_dev $RANDOM $DYN_RAN` #normal distibution

# # # change the sample from log scale back to natural scale
# # set sd_sample=`echo "" | awk '{print ((10^('$sd_sample_log10')))}'`
# # echo "a sd sample= "$sd_sample