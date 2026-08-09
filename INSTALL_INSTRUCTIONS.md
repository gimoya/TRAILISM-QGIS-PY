# How to Make Script Visible in QGIS Processing Toolbox

## Quick Steps:

1. **Find your QGIS scripts folder:**
   - Open QGIS
   - Go to: **Settings → Options → Processing → Scripts**
   - Note the folder path shown (e.g., `C:\Users\YourName\AppData\Roaming\QGIS\QGIS3\profiles\default\processing\scripts`)

2. **Copy the script:**
   - Copy `TRAILISM__Add_Trail_Shapes.py` to that folder

3. **Reload scripts:**
   - In QGIS: **Processing → Toolbox**
   - Right-click on "Scripts" → **Refresh Scripts**
   - OR simply restart QGIS

4. **Find the algorithm:**
   - In Processing Toolbox, look under: **Scripts → TRAILISM → TRAILISM: Add Trail Shapes**
   - Or search for "Add Trail Shapes" in the toolbox search bar

## Alternative: Check if script is already in scripts folder

If the script is already there but not showing:
- Make sure the file name is exactly `TRAILISM__Add_Trail_Shapes.py`
- Check for syntax errors (the script should compile without errors)
- Try restarting QGIS completely

## View Debug Output

When you run the algorithm, the debug output will appear in:
- **Processing Log** panel (View → Panels → Log Messages Panel)
- Or at the bottom of the Processing Toolbox window in the "Log" tab

