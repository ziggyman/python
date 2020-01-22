#!/usr/bin/env bash
# vi: ts=3:
# Adapted from shocdatacopy.sh
#Version 0.9
# Version 2.2
# 03/02/2016

# ----- Global variables

DATE_ENTERED=0
VERSION=0.9

PATH=$PATH:/usr/local/bin
LOGPATH=/var/log/datatransfer/
INST_NAME=$(hostname)
SHORTNAME=spupnic
DATA_PREFIX=/home/ccd/data

if [ -z $USER ] 
then
	[ -z $LOGNAME ] && exit 1 || USER=$LOGNAME
fi

# ----- Global variables

TELESCOPE='74in'

# ---- Functions used in this script ----

function next_in_seq () {
	# Function to work out what the latest log filename is and add 1 to 
	# that for the day. This will ensure that in time, the log files are 
	# overwritten but they'll be there for the astronomer to check
	# if they want

	FILENAME=$1
	LATEST_FILE=$(ls -lrt ${FILENAME}* 2>/dev/null|tail -1|cut -d'/' -f4|cut -d'-' -f2)
	[ "${LATEST_FILE}x" == "x" ] && NEXT_SEQ=1 || NEXT_SEQ=$((${LATEST_FILE} + 1))
	export NEXT_SEQ
	return 
}

# Work out whether we have changed year or day
function ymd () {
	# If there is NO directory for today, then we must be transferring
	# yesterday's data.
	# If there's a change in year, we need to go back 1 year and transfer
	# data from the last day of last year
#	if ! [ -d ${DATA_PREFIX}/${TELESCOPE}/${SHORTNAME}/${YYYY}/${MMDD} ]
    # TODO: change this to standard path
	if ! [ -d ${DATA_PREFIX}/${YYYY}${MMDD} ]
	then
	# If the directory doesn't exist, then we must be into the next morning
	# and we should not try to transfer the data from today but rather from yesterday
	# Let's check?
		TM=`date +'%H%M'`
		if [ ${TM} -gt 0000 ] && [ ${TM} -lt 1100 ]
		then
		# It's between midnight and 11am, which means any time this script is being
		# run now, it's probably wanting to copy last night's data
		# If this is not what you want, then run the script with a '-d' option
			MMDD=$(date --date='Yesterday' +'%m%d')
			YYYY=$(date --date='Yesterday' +'%Y')
		fi
	fi
	return
}

# ---- END Functions used in this script ----

# Get options from the command line if there are any

while getopts hd: opts
do
	case $opts in
		d|D)	DT=${OPTARG}
				YYYYMMDD=$(echo $DT|tr '[A-Z]' '[a-z]')
				DATE_ENTERED=1
				
				if [ ${YYYYMMDD} = 'yesterday' ] 
				then
					DTYYYY=$(date --date='Yesterday' +'%Y')
					DTMMDD=$(date --date='Yesterday' +'%m%d')
#					if ! [ -d ${DATA_PREFIX}/${TELESCOPE}/${SHORTNAME}/${DTYYYY}/${DTMMDD} ]
					if ! [ -d ${DATA_PREFIX}/${DTYYYY}${DTMMDD} ]
					then
#						echo "Sorry. No directory for yesterday (${SHORTNAME}/$DTYYYY${DTMMDD}) to copy"
						echo "Sorry. No directory for yesterday (${DATA_PREFIX}/${DTYYYY}${DTMMDD}) to copy"
						exit 1
					fi
				elif [ ${YYYYMMDD} = 'today' ]
				then
					DTYYYY=$(date --date='Today' +'%Y')
					DTMMDD=$(date --date='Today' +'%m%d')
#					if ! [ -d ${DATA_PREFIX}/${TELESCOPE}/${SHORTNAME}/${DTYYYY}/${DTMMDD} ]
					if ! [ -d ${DATA_PREFIX}/${DTYYYY}${DTMMDD} ]
					then
						echo "Sorry. No directory for today (${DATA_PREFIX}/$DTYYYY${DTMMDD}) to copy"
						exit 1
					fi
				else
					echo "Date entered was ${YYYYMMDD}"
					DTYYYY=${YYYYMMDD:0:4}
					DTMMDD=${YYYYMMDD:4}
#					if ! [ -d ${DATA_PREFIX}/${TELESCOPE}/${SHORTNAME}/${DTYYYY}/${DTMMDD} ]
					if ! [ -d ${DATA_PREFIX}/${DTYYYY}${DTMMDD} ]
					then
						echo "Sorry. No directory (${SHORTNAME}/$DTYYYY${DTMMDD}) to transfer"
						exit 1
					fi
			fi

				;;
		h|H) echo "Usage: $0 [ -d YYYYMMDD|today|yesterday ] [-h]"
				echo "Script version: ${VERSION}"
				exit 0
				;;
	esac
done

# ---- Build the logfile ----
FILE="${TELESCOPE}DataCopyLog"

MMDD=`date +'%m%d'`
next_in_seq ${LOGPATH}/${FILE}${MMDD}

LOG="${TELESCOPE}DataCopyLog${MMDD}-${NEXT_SEQ}"
if [ -e ${LOGPATH}/${LOG} ]
then
	echo "Found that this logfile (${LOGPATH}/${LOG} is already there"
	echo "Creating alternative by appending time to log file"
	HHMM=$(date +'%H%M')
	LOG="${LOG}.${HHMM}"
fi
# ---- End of building the logfile ----

echo "Use an editor/less to view the logfile ${LOGPATH}/${LOG} for output of this data copy"

exec 2>${LOGPATH}/${LOG}
exec 1>&2

echo -e "User: $USER"

echo -e "Output log of this data copy in ${LOGPATH}/${LOG}\n\t"

#echo "Trying to work out whether I'm in Sutherland or Cape Town"

SL=astro.suth.saao.ac.za
CT=astro.cape.saao.ac.za

DEF_ROUTE=$(ip route list|grep default|cut -d' ' -f3)

if [ "${DEF_ROUTE}x" == '10.2.1.1x' ]
then
	DEST_HOST=${SL} 
	echo "I'm in Sutherland, copying data to ${DEST_HOST}"
else 
	DEST_HOST=${CT}
	echo "I'm in Cape Town, copying data to ${DEST_HOST}"
fi

HHMM=`date +'%H:%M'`

if [ $DATE_ENTERED -eq 1 ]
then
	echo -e "** Copying data (as per your request), ${YYYYMMDD} @ $HHMM\n"
	YYYY=${DTYYYY}
	MMDD=${DTMMDD}
else
	YYYY=$(date +'%Y')
	MMDD=$(date +'%m%d')
	# Call the date function to work out whether we're trying to transfer
	# yesterday's data this morning
	ymd
fi


[ $DATE_ENTERED -ne 1 ] && echo -e "** Copying data ($YYYY), $MMDD @ $HHMM\n"
#echo "Would call rsync -uavz ${DATA_PREFIX}/$TELESCOPE/${SHORTNAME}/${YYYY}/${MMDD} ${DEST_HOST}::saaotelescopedata/$TELESCOPE/$SHORTNAME/${YYYY}"
echo "Calling rsync -uavz ${DATA_PREFIX}/${YYYY}${MMDD} ${DEST_HOST}::saaotelescopedata/$TELESCOPE/$SHORTNAME/${YYYY}/${MMDD}"
#time rsync -uavz ${DATA_PREFIX}/$TELESCOPE/${SHORTNAME}/${YYYY}/${MMDD} ${DEST_HOST}::saaotelescopedata/$TELESCOPE/$SHORTNAME/${YYYY}
time rsync -uavz ${DATA_PREFIX}/${YYYY}${MMDD} ${DEST_HOST}::saaotelescopedata/$TELESCOPE/$SHORTNAME/${YYYY}/${MMDD}

[ $? -ne 0 ] && echo "There were no files to transfer for ${YYYY}${MMDD}" 

echo -e "\n** Done"
