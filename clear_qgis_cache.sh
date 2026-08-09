#!/bin/bash
# Clear QGIS cache script

APPDATA_PATH="/c/Users/kayci_76pke/AppData/Roaming/QGIS/QGIS3/cache"
LOCALAPPDATA_PATH="/c/Users/kayci_76pke/AppData/Local/QGIS/QGIS3/cache"

echo "Clearing QGIS cache..."

if [ -d "$APPDATA_PATH" ]; then
    echo "Found cache in Roaming: $APPDATA_PATH"
    rm -rf "$APPDATA_PATH"/* 2>/dev/null
    if [ $? -eq 0 ]; then
        echo "Deleted Roaming cache"
    else
        echo "Roaming cache already empty or error deleting"
    fi
else
    echo "Roaming cache directory not found"
fi

if [ -d "$LOCALAPPDATA_PATH" ]; then
    echo "Found cache in Local: $LOCALAPPDATA_PATH"
    rm -rf "$LOCALAPPDATA_PATH"/* 2>/dev/null
    if [ $? -eq 0 ]; then
        echo "Deleted Local cache"
    else
        echo "Local cache already empty or error deleting"
    fi
else
    echo "Local cache directory not found"
fi

echo "Done! Please restart QGIS."

