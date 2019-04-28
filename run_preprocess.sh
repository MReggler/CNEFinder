#!/bin/bash

# Rough source:
#https://gist.github.com/clcollins/68cd6b5a67a6124cafbf


# Set to true to just see what command
# would be run
DRY_RUN=false

IMAGE='preprocess'
# Local volume to be mapped into the container any time you run it
# usually with config files or whatnot
MYCONFDIR=""

# An array of envvars
ENVVARS=('APPENV docker')

# prevent creation of subshell by `not` using a pipe
# https://unix.stackexchange.com/questions/99423
file=".env"
while read line; do
    ENVVARS+=($line)
done <<<"$(sed -e '/^\s*$/ d' -e '/^#/ d' $file)"

# Array of ports like (80:80 443:443 3000 8080)
# Can be mapped or unmapped
PORTS=()

# Array of volumes
VOLUMES=("/cnefinder/input:/input" "/cnefinder/output:/output")
ENTRYPOINT="" #"/bin/bash"
CMD=""
RESTART=""
DAEMON=false

# The Docker command to use.
# Could be different if including --tlsverify -H "hostname:hostport", etc
DOCKERCMD="docker"
NAME=${IMAGE}

declare -a ENVVAR_STRING
for envvar in ${ENVVARS[@]} ; do
  ENVVAR_STRING+=("-e ${envvar}")
done

declare -a PORT_STRING
for port in ${PORTS[@]} ; do
  PORT_STRING+=("-p ${port}")
done

declare -a VOLUME_STRING
for volume in ${VOLUMES[@]} ; do
  VOLUME_STRING+=("-v ${volume}")
done

if [[ ! -z $NAME ]] ; then
  NAME_STRING="--name ${NAME}"
fi

if [[ ! -z $RESTART ]] ; then
  RESTART_STRING="--restart ${RESTART}"
fi

if [[ ! -z $ENTRYPOINT ]] ; then
  ENTRYPOINT_STRING="--entrypoint ${ENTRYPOINT}"
fi

if [[ ! -z $CMD ]] ; then
  CMD_STRING="${CMD}"
fi

if $DAEMON ; then
  DAEMON_STRING='-d'
else
  DAEMON_STRING='-it'
fi

if $DOCKERCMD ps -a | awk "/${NAME}/ {print $NF}" | grep $NAME &>/dev/null ; then
  $DOCKERCMD stop $NAME 1>/dev/null \
  && $DOCKERCMD rm $NAME 1>/dev/null
fi

OPTS="${ENVVAR_STRING[@]} ${PORT_STRING[@]} ${VOLUME_STRING[@]} $NAME_STRING $RESTART_STRING $ENTRYPOINT_STRING"

THE_COMMAND="$DOCKERCMD run $OPTS $DAEMON_STRING $IMAGE $CMD_STRING"
if $DRY_RUN ; then
  echo "$THE_COMMAND"
  exit 0
else
  echo "$THE_COMMAND"
  exec $THE_COMMAND
fi
