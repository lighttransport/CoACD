// SPDX-License-Identifier: Apache-2.0
// Thin ES6 wrapper around the CoACD WASM module.
//
//   import { loadCoACD } from './CoACDLoader.js';
//   const coacd = await loadCoACD();
//   const hulls = coacd.decompose(vertices, indices, { threshold: 0.1 });

let _modulePromise = null;

export async function loadCoACD(opts = {}) {
    if (!_modulePromise) {
        const url = opts.scriptURL ?? new URL('./src/coacd.js', import.meta.url).href;
        const mod = await import(/* @vite-ignore */ url);
        _modulePromise = mod.default(opts.moduleArg ?? {});
    }
    const Module = await _modulePromise;

    return {
        /**
         * @param {Float32Array|Float64Array} vertices  flat xyz
         * @param {Int32Array|Uint32Array} indices      flat triangle indices
         * @param {object} [options]
         * @returns {Array<{vertices: Float64Array, indices: Int32Array}>}
         */
        decompose(vertices, indices, options = {}) {
            const v64 = vertices instanceof Float64Array
                ? vertices : new Float64Array(vertices);
            const i32 = indices instanceof Int32Array
                ? indices : new Int32Array(indices);
            return Module.decompose({
                vertices: v64,
                indices: i32,
                options,
            });
        },
        setLogLevel(level) {
            Module.set_log_level(level);
        },
        _raw: Module,
    };
}
