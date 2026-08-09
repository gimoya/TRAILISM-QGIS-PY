# Auto-Stage and Auto-Commit Setup

## Scripts overview

| Script | Use when |
|--------|----------|
| **auto_stage_commit.sh** | Watcher: auto-stage and commit on a configurable interval (pattern + interval in seconds). |
| **one_click_push.sh** | One-click: stage all, commit, push to `main`; then waits for Enter. |

## Quick Start

**Run the auto-stage-commit watcher:**
```bash
./auto_stage_commit.sh
```

This will:
- Auto-stage any changes to `*.py` files
- Auto-commit every **300 seconds (5 minutes)** if there are staged changes (pass a second argument to override, e.g. `./auto_stage_commit.sh "*.py" 60` for 1 minute)
- Create a complete git history log

## Options

**Watch specific files:**
```bash
./auto_stage_commit.sh "TRAILISM__Add_Trail_Shapes.py" 30
# Watches only that file, commits every 30 seconds
```

**Watch all Python files with custom interval:**
```bash
./auto_stage_commit.sh "*.py" 120
# Commits every 2 minutes
```

## What Gets Committed

- All changes to matching files are automatically staged
- Commits are created with timestamps: `Auto-commit #N: YYYY-MM-DD HH:MM:SS`
- Full git history is maintained - you can always compare versions

## View History

```bash
# See all auto-commits
git log --oneline --grep="Auto-commit"

# Compare two commits
git diff <commit1> <commit2> -- TRAILISM__Add_Trail_Shapes.py

# See what changed in a file over time
git log --oneline -- TRAILISM__Add_Trail_Shapes.py
```

## Stop the Watcher

Press `Ctrl+C` to stop the auto-commit watcher.

## Notes

- The watcher runs in the foreground (you'll see output)
- Commits only happen if there are actual changes
- Each commit creates a snapshot you can compare later
- Safe to run - won't overwrite uncommitted work
