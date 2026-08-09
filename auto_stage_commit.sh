#!/bin/bash
# Auto-stage and auto-commit script
# Watches for file changes, auto-stages them, and commits to create history
# Usage: ./auto_stage_commit.sh [file_pattern] [commit_interval_seconds]

cd "$(dirname "$0")"

PATTERN="${1:-*.py}"
INTERVAL="${2:-300}"  # Default: commit every 300 seconds (5 minutes) if changes detected

echo "Auto-staging and committing changes to: $PATTERN"
echo "Commit interval: $INTERVAL seconds"
echo "Press Ctrl+C to stop"
echo ""

LAST_COMMIT=$(date +%s)
COMMIT_COUNT=0

while true; do
    sleep 5  # Check every 5 seconds
    
    # Check if files matching pattern have changed
    CHANGED_FILES=$(git diff --name-only "$PATTERN" 2>/dev/null)
    STAGED_FILES=$(git diff --cached --name-only "$PATTERN" 2>/dev/null)
    
    if [ -z "$CHANGED_FILES" ] && [ -z "$STAGED_FILES" ]; then
        continue  # No changes
    fi
    
    # Auto-stage all changes matching pattern
    if [ -n "$CHANGED_FILES" ]; then
        echo "[$(date +'%H:%M:%S')] Auto-staging: $CHANGED_FILES"
        git add "$PATTERN" 2>/dev/null
    fi
    
    # Check if we have staged changes
    if git diff --cached --quiet; then
        continue  # Nothing staged
    fi
    
    CURRENT_TIME=$(date +%s)
    TIME_SINCE_COMMIT=$((CURRENT_TIME - LAST_COMMIT))
    
    # Commit if enough time has passed
    if [ $TIME_SINCE_COMMIT -ge $INTERVAL ]; then
        COMMIT_COUNT=$((COMMIT_COUNT + 1))
        TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')
        COMMIT_MSG="Auto-commit #$COMMIT_COUNT: $TIMESTAMP"
        
        git commit -m "$COMMIT_MSG" > /dev/null 2>&1
        
        if [ $? -eq 0 ]; then
            LAST_COMMIT=$(date +%s)
            STAGED_COUNT=$(git diff --cached --name-only | wc -l)
            echo "[$(date +'%H:%M:%S')] ✓ Committed: $COMMIT_MSG"
            echo "   Files: $(git log -1 --name-only --pretty=format:'' | head -5 | tr '\n' ' ')"
            echo ""
        fi
    else
        REMAINING=$((INTERVAL - TIME_SINCE_COMMIT))
        echo "[$(date +'%H:%M:%S')] Staged changes waiting... (commit in ${REMAINING}s)"
    fi
done
