const { contextBridge, ipcRenderer } = require('electron');
contextBridge.exposeInMainWorld('relibDesktop', {
  version: () => ipcRenderer.invoke('version'),
  check: () => ipcRenderer.invoke('check'),
  download: () => ipcRenderer.invoke('download'),
  install: () => ipcRenderer.invoke('install'),
  savePdf: (name, data) => ipcRenderer.invoke('savePdf', { name, data }),
  onUpdate: cb => ipcRenderer.on('upd', (e, m) => cb(m))
});
