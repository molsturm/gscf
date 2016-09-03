#!/bin/sh
## ---------------------------------------------------------------------
## Put licence text here
## ---------------------------------------------------------------------
#
# Settings
#
# The git repo to checkout
FROM="https://github.com/linalgwrap/linalgwrap"

# Branch to checkout (empty for default)
BRANCH=""

# Select this branch
BRANCH="master-krims"

# Folder to check it out to
WHAT="linalgwrap"

# Interval: How often to update:
INTERVAL="1 hour"

# File to use in order to test a successful checkout
CHECKFILE="CMakeLists.txt"

if [ ! -f "$PWD/get.lib.sh" ]; then
	echo "PWD needs to be the location of the get.lib.sh file."
	exit 1
fi

. "$PWD/get.lib.sh" || exit 1
update_repo "$FROM" "$WHAT" "$CHECKFILE" "$INTERVAL" "$BRANCH"
exit $?
