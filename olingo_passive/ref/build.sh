#!/bin/bash

users=""
threshold=""

for arg in "$@"; do
	case $arg in
		--users=*)
			users="${arg#*=}"
			shift
			;;
		--threshold=*)
			threshold="${arg#*=}"
			shift
			;;
		*)
			echo "Unknown option $arg"
			exit 1
			;;
	esac
done

if [ -z "$users" ] || [ -z "$threshold" ]; then
    echo "Usage: $0 --users=<number> --threshold=<number>"
    exit 1
fi

echo "USERS: $users"
echo "THRESHOLD: $threshold"

if ! [[ "$users" =~ ^[0-9]+$ ]] || ! [[ "$threshold" =~ ^[0-9]+$ ]]; then
    echo "Error: Both users and threshold must be numbers"
    exit 1
fi

make bench USERS=$users THRESHOLD=$threshold DKG=1
make bench USERS=$users THRESHOLD=$threshold PREENC=1
make bench USERS=$users THRESHOLD=$threshold ENCRYPT=1
make bench USERS=$users THRESHOLD=$threshold AS2=1
make ghkss USERS=$users THRESHOLD=$threshold
