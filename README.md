## Fast Fourier That For You

This is an in-place, radix 2 implementation of the Cooley-Tukey FFT algorithm.

## Usage

```lua
local fft = require("fft")

local SIZE = 32

local sine = {r = table.create(SIZE, 0), c = table.create(SIZE, 0) }
for i = 1, SIZE do
  sine.r[i] = math.sin(4*math.pi*i / SIZE)
end

local X = fft(sine)

local out = {}
for i = 1, SIZE do
  out[i] = string.format("(%.2f,%.2f)", X.r[i], X.c[i]) .. (if i % 8 == 0 then "\n" else "")
end
print("{" .. table.concat(out, ",") .. "}")

local x = fft(X, true, false, -1)

print("R: recovered input; I: initial input")
for i = 1, SIZE do
  out[i] = string.format("R:%.2f, I:%.2f", x.r[i], sine.r[i]) .. (if i % 2 == 0 then "\n" else "")
end
print(table.concat(out, "; "))

--[=[
{(-0.00,0.00),(-0.00,-16.00),(0.00,-0.00),(-0.00,-0.00),(-0.00,-0.00),(-0.00,-0.00),(0.00,-0.00),(0.00,-0.00)
,(0.00,0.00),(0.00,-0.00),(-0.00,0.00),(0.00,-0.00),(0.00,0.00),(-0.00,0.00),(-0.00,0.00),(-0.00,0.00)
,(-0.00,-0.00),(-0.00,-0.00),(0.00,-0.00),(0.00,0.00),(-0.00,-0.00),(0.00,0.00),(0.00,-0.00),(0.00,0.00)
,(0.00,0.00),(-0.00,0.00),(-0.00,0.00),(-0.00,0.00),(0.00,0.00),(-0.00,16.00),(-0.00,-0.00),(-0.00,0.00)
}
R: recovered input; I: initial input
R:0.38, I:0.38; R:0.71, I:0.71
; R:0.92, I:0.92; R:1.00, I:1.00
; R:0.92, I:0.92; R:0.71, I:0.71
; R:0.38, I:0.38; R:0.00, I:0.00
; R:-0.38, I:-0.38; R:-0.71, I:-0.71
; R:-0.92, I:-0.92; R:-1.00, I:-1.00
; R:-0.92, I:-0.92; R:-0.71, I:-0.71
; R:-0.38, I:-0.38; R:-0.00, I:-0.00
; R:0.38, I:0.38; R:0.71, I:0.71
; R:0.92, I:0.92; R:1.00, I:1.00
; R:0.92, I:0.92; R:0.71, I:0.71
; R:0.38, I:0.38; R:0.00, I:0.00
; R:-0.38, I:-0.38; R:-0.71, I:-0.71
; R:-0.92, I:-0.92; R:-1.00, I:-1.00
; R:-0.92, I:-0.92; R:-0.71, I:-0.71
; R:-0.38, I:-0.38; R:-0.00, I:-0.00
]=]
```

## Documentation

Input should be of the form: `type zvec = { r: { number }, c: { number } }`

Input & output vector size will be zero-padded to the next power of two wrt real input vector size.

There are no dependencies.

### Arguments to `fft`

- `zvec x`: The discrete sample domain.

- `boolean? invert`: Inverse DFT is computed by passing `true`.

- `boolean? keep_revert`: The FFT algorithm produces the DFT in reverse, and then re-reverses it later. You can avoid that small overhead by passing `true`. The output zvec will be backwards.

- `number? scale_pow`: This parameter controls the exponent of a scalar that is multiplied on the output. The base is the output size. Passing `-1` normalizes the output with `1/N`, and `-0.5` with `1/sqrt(N)`.


