-- @densenet
-- SPDX-License-Identifier: Apache-2.0

local TAU = math.pi*2

local lz, l = bit32.countlz, bit32.lshift

local cos, sin = math.cos, math.sin

local make, copy = table.create, table.clone

local __rev = {}

export type zvec = { r: {number}, c: {number} }

local function fft(x: zvec, invert: boolean?, keep_revert: boolean?, scale_pow: number?): zvec
  local N_ORIG = #x.r
  if N_ORIG <= 1 then
    return {r=copy(x.r),c=copy(x.c)}
  end
  local K = (32 - lz(N_ORIG-1))
  local N = l(1, K)
  
  local X_r: {number}, X_c: {number} = make(N), make(N)
  
  -- generate the Kth bit-reversal permutation
  local rev = __rev[K]
  if rev == nil then
    rev = make(N, 0)
    __rev[K] = rev
    do
      local n0 = 0
      local n = 1
      for k = 1, K do
        n0 = n
        n *= 2
        for i = 1, n0 do
          rev[i] = rev[i]*2
        end
        for i = n0+1, n do
          rev[i] = rev[i-n0]+1
        end
      end
    end
  end
  
  -- permuted copy
  for k = 1, N do
    X_r[rev[k]+1], X_c[rev[k]+1] = x.r[k], x.c[k]
  end  
  local si = if not invert then -1 else 1
  
  -- cooley tukey
  do
    local m = 1
    local wm_r, wm_c, w_r, w_c
    local t_r, t_c, u_r, u_c
    local i0, i1, X_r0, X_c0
    for s = 1, K do
      m *= 2
      wm_r, wm_c = cos(TAU/m), si*sin(TAU/m)
      for k = 0, N-1, m do
        w_r, w_c = 1, 0
        for j = 0, m/2-1 do
          i0, i1 = 1+k+j, 1+k+j+m/2
          X_r0, X_c0 = X_r[i0], X_c[i0]
          u_r, u_c = X_r[i1], X_c[i1]
        
          t_r, t_c = w_r * X_r0 - w_c * X_c0, w_r * X_c0 + w_c * X_r0
        
          X_r[i0], X_c[i0] = u_r + t_r, u_c + t_c
          X_r[i1], X_c[i1] = u_r - t_r, u_c - t_c
        
          w_r, w_c = w_r * wm_r - w_c * wm_c, w_r * wm_c + w_c * wm_r
        end
      end
    end
  end
  
  if not keep_revert then
    for i = 1, N/2 do
      X_r[i], X_c[i], X_r[1+N - i], X_c[1+N - i] = X_r[1+N - i], X_c[1+N - i], X_r[i], X_c[i]
    end
  end
  
  if scale_pow then
    local scalar = if scale_pow == -1 then 1/N elseif scale_pow == -0.5 then 1/math.sqrt(N) else math.pow(N, scale_pow)
    for i = 1, N, 2 do
      X_r[i], X_c[i], X_r[i+1], X_c[i+1] = scalar*X_r[i], scalar*X_c[i], scalar*X_r[i+1], scalar*X_c[i+1]
    end
  end
  
  return {r=X_r,c=X_c}
end

local __DEBUG = nil

local function __dbg(env)
  local function __dbg_print(x)
    local buf = {}
    
    for i, v in x.r do
      buf[#buf+1] = string.format("[%d]=(%.3f, %.3f)", i, x.r[i], x.c[i])
    end
    
    print(string.format("{ %s }", table.concat(buf, ", ")))
  end

  local function __dbg_gen(N)
    local x_r, x_c = table.create(N), table.create(N)
  
    for i = 1, N do
      x_r[i], x_c[i] = math.random(0, 255), math.random(0, 255)
    end
  
    return {r = x_r, c = x_c}
  end
  
  local function __dbg_gen_sine(N, freq)
    local x_r, x_c = table.create(N), table.create(N, 0)
    
    for i = 1, N do
      x_r[i] = math.cos(freq*math.pi*(i/N))
--      x_c[i] = math.sin(freq*math.pi*(i/N))
    end
    
    return {r=x_r,c=x_c}
  end
  
  if __DEBUG == 1 then
    local x = __dbg_gen(15)
    __dbg_print(fft(x))
  end
  
  if __DEBUG == 2 then
    local x = __dbg_gen_sine(env.N or 64, env.A or 1)
    __dbg_print(fft(x))
  end
  
  if __DEBUG >= 3 and __DEBUG < 5 then
    local tenx = require(".\\vendor\\luau-fft\\fourier")
    local tenxc = require(".\\vendor\\luau-fft\\cmath")
  
    if __DEBUG == 3 then
      local x = __dbg_gen_sine(env.N or 64, env.A or 1)
    
      local X_tenx = tenx.fft((function()
        local x2cmath = {}
      
        for i, v in x.r do
          x2cmath[i] = tenxc.new(v, x.c[i])
        end
      
        return x2cmath
      end)())
    
      local X_loc = (function()
        local X = {r={},c={}}
      
        for i, c in X_tenx do
          X.r[i], X.c[i] = c.Real, c.Imag
        end
      
        return X
      end)()
    
      __dbg_print(X_loc)
    end
  
    if __DEBUG == 4 then
      local x = __dbg_gen_sine(64, 2)
      local x_tenx = (function()
        local x2cmath = {}
      
        for i, v in x.r do
          x2cmath[i] = tenxc.new(v, x.c[i])
        end
      
        return x2cmath
      end)()
      local fft_tenx = tenx.fft
      
      local stamp = os.clock()
      
      local rep = 1e4
      
      for i = 1, rep do
        fft_tenx(x_tenx)
      end
      
      print(string.format("tenx: %.3fms", 1e3*(os.clock()-stamp)))
      
      local stamp = os.clock()
      
      for i = 1, rep do
        fft(x)
      end
      
      print(string.format("local: %.3fms", 1e3*(os.clock()-stamp)))
    end
  end
end

--__DEBUG = 3
--__dbg {N=32,A=2}

if not __DEBUG then
  return fft
end