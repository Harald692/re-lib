// re.lib – Windows-App (Electron)
const { app, BrowserWindow, ipcMain, dialog, shell, Menu } = require('electron');
const path = require('path');
const fs = require('fs');
const { autoUpdater } = require('electron-updater');

autoUpdater.autoDownload = false;          // erst nach Klick auf "Jetzt aktualisieren"
autoUpdater.autoInstallOnAppQuit = true;
let win = null;

if (!app.requestSingleInstanceLock()) app.quit();
app.on('second-instance', () => { if (win) { if (win.isMinimized()) win.restore(); win.focus(); } });

function createWindow() {
  win = new BrowserWindow({
    width: 1280, height: 900, minWidth: 380, minHeight: 500,
    title: 're.lib', backgroundColor: '#173a6b',
    icon: path.join(__dirname, 'www', 'icon-512.png'),
    webPreferences: { preload: path.join(__dirname, 'preload.js'), contextIsolation: true, nodeIntegration: false }
  });
  Menu.setApplicationMenu(null);
  win.loadFile(path.join(__dirname, 'www', 'index.html'));
  // Links nach aussen im Standardbrowser öffnen
  win.webContents.setWindowOpenHandler(({ url }) => { if (/^https?:/.test(url)) shell.openExternal(url); return { action: 'deny' }; });
}
app.whenReady().then(createWindow);
app.on('window-all-closed', () => app.quit());

// ---- Updates über GitHub-Releases ----
const send = (s, d) => { if (win) win.webContents.send('upd', { s, d }); };
autoUpdater.on('update-available', i => send('available', i.version));
autoUpdater.on('update-not-available', () => send('none'));
autoUpdater.on('download-progress', p => send('progress', Math.round(p.percent)));
autoUpdater.on('update-downloaded', () => send('ready'));
autoUpdater.on('error', e => send('error', String((e && e.message) || e)));

ipcMain.handle('version', () => app.getVersion());
ipcMain.handle('check', async () => {
  if (!app.isPackaged) return { dev: true };
  try { await autoUpdater.checkForUpdates(); } catch (e) { send('error', String(e.message || e)); }
  return {};
});
ipcMain.handle('download', () => autoUpdater.downloadUpdate());
ipcMain.handle('install', () => { setImmediate(() => autoUpdater.quitAndInstall(false, true)); });

// ---- PDF speichern und öffnen (zum Drucken) ----
ipcMain.handle('savePdf', async (e, { name, data }) => {
  const r = await dialog.showSaveDialog(win, {
    title: 'PDF speichern', defaultPath: path.join(app.getPath('documents'), name),
    filters: [{ name: 'PDF-Dokument', extensions: ['pdf'] }]
  });
  if (r.canceled || !r.filePath) return null;
  fs.writeFileSync(r.filePath, Buffer.from(data, 'base64'));
  shell.openPath(r.filePath);
  return r.filePath;
});
