#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <file_extension>"
    exit 1
fi

DIRECTORY=$1
EXTENSION=$2

# Find all files with the given extension recursively and trim trailing whitespace
find "$DIRECTORY" -type f -name "*.$EXTENSION" | while read -r file; do
    sed -i '' -e 's/[ \t]*$//' "$file"
done

echo "Trailing whitespace trimmed for all *.$EXTENSION files in $DIRECTORY"
