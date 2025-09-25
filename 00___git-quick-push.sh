#!/bin/bash

# Usage:
#   ./git-quick-push.sh
#   ./git-quick-push.sh "Your commit message here"

MSG="${*:-latest updates}"

echo "Staging changes..."
git add -A

echo "Committing (\"$MSG\")..."
if git commit -m "$MSG" 2>/dev/null; then
    echo "Committed successfully."
else
    echo "Nothing to commit or commit skipped; continuing."
fi

echo "Pushing..."
git push -u origin master

read -p "Press [Enter] to continue..."
