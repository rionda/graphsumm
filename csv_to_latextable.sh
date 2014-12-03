#! /bin/sh

if [ $# -lt 1 ]; then
  echo "USAGE: $0 csvfile" > /dev/stderr
  exit 1
fi

if [ ! -r $1 ]; then
  echo "ERROR: Input file $1 not readable" > /dev/stderr
  exit 1
fi

sed -e 's/,/ \& /g' -e 's/$/ \\\\/' $1

