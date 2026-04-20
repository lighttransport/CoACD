# CoACD WASM/Emscripten build

Single-threaded browser/Node WebAssembly build of CoACD with embind bindings.

## Build

```bash
# requires emcmake on PATH
./bootstrap-linux.sh
# or:  BUILD_TYPE=Debug ./bootstrap-linux.sh
```

Artifacts land in `js/src/coacd.{js,wasm}`.

## Usage (browser)

```js
import { loadCoACD } from './js/CoACDLoader.js';
const coacd = await loadCoACD();
const hulls = coacd.decompose(verticesF64, indicesI32, { threshold: 0.05 });
```

Serve the `web/` directory (e.g. `python3 -m http.server`) and open
`demo/index.html`.

## Node smoke test

```bash
node --experimental-wasm-bigint -e "
  import('./js/src/coacd.js').then(async m => {
    const M = await m.default();
    const v = new Float64Array([0,0,0, 1,0,0, 0,1,0, 0,0,1]);
    const f = new Int32Array([0,2,1, 0,1,3, 0,3,2, 1,2,3]);
    const r = M.decompose({ vertices: v, indices: f, options: { threshold: 0.2, preprocess: 'off' } });
    console.log('hulls:', r.length);
  });
"
```

## Limitations

- Single-threaded (no pthreads, no OpenMP).
- OpenVDB-based preprocessing is **disabled**; pass `preprocess: "off"`
  (the default). `"on"` / `"auto"` will not work.
- CDT triangulation is disabled; the built-in ear-clipping path is used.
- Inputs/outputs are typed arrays, not `coacd::Mesh`.
