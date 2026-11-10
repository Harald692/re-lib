const C='relib-__BUILD__';
const F=['./','./index.html','./manifest.webmanifest','./icon-192.png','./icon-512.png','./icon-maskable-512.png','./jspdf.umd.min.js','./pdffont.js'];
self.addEventListener('install',e=>{e.waitUntil(caches.open(C).then(c=>c.addAll(F)));self.skipWaiting()});
self.addEventListener('activate',e=>{e.waitUntil(caches.keys().then(k=>Promise.all(k.filter(x=>x!==C).map(x=>caches.delete(x)))));self.clients.claim()});
self.addEventListener('fetch',e=>{
  if(e.request.url.includes('version.json'))return; // immer frisch vom Server
  e.respondWith(caches.match(e.request).then(r=>r||fetch(e.request)));
});
