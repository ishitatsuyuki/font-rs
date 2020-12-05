import * as wasm from "font-tuner";

async function update() {
    if(document.getElementById('filepicker').files.length === 0) return
    const file = document.getElementById('filepicker').files[0]
    const text = document.getElementById('textinput').value
    const buf = wasm.render(new Uint8Array(await file.arrayBuffer()), text)
    document.getElementById('canvas').getContext('2d').putImageData(new ImageData(new Uint8ClampedArray(buf), 512), 0, 0)
}

document.getElementById('filepicker').addEventListener('change', update, false)
document.getElementById('textinput').addEventListener('input', update)