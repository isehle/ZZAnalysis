#!/bin/tcsh -f

# parameter: XX = number of the bigbirdXX schedd; to look at loads in the different schedds run:
# condor_status -schedd

# Make the grid proxy available in ~, if existing and valid
set proxy_valid=`voms-proxy-info --timeleft`
if ($proxy_valid > 10 ) then
   echo "GRID proxy found, validity: $proxy_valid s"
   if ($?X509_USER_PROXY) then
    if ($X509_USER_PROXY != ~/x509up_u${uid}) cp $X509_USER_PROXY ~/x509up_u${uid}
   else if (-e /tmp/x509up_u${uid} ) then
     cp /tmp/x509up_u${uid} ~
   endif
else # Last attempt: Try to see if a valid proxy already exists in ~
   if ( -e ~/x509up_u${uid} ) then
      setenv X509_USER_PROXY ~/x509up_u${uid}
      set proxy_valid=`voms-proxy-info --timeleft`
      if ($proxy_valid > 10 ) then
         echo "GRID proxy found in ~, validity: $proxy_valid s"
      endif
   endif

   if ($proxy_valid < 10 ) then
      echo "Error: no valid GRID proxy found."
      exit 1
   endif
endif

set JOBNAME=`basename $PWD`

set queue=' -queue directory in'

foreach x (*Chunk*)
 set nonomatch
 set logFile = ( ${x}/log/*.log )
 if ( -e $logFile[1] ) then
    echo "\n${x}: job already submitted. If you want to resubmit, ensure all jobs are finished and run cleanup.csh.\nAborting."
    exit 1
 endif
 unset nonomatch


 set queue="$queue $x"
end

if (! ($1 == "") ) then
   setenv _CONDOR_SCHEDD_HOST bigbird$1.cern.ch
   setenv _CONDOR_CREDD_HOST bigbird$1.cern.ch
endif

condor_submit condor.sub $queue
