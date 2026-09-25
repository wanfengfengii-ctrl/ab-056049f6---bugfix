// 纯计算几何工具：二维点、有符号距离、凸包、点与简单多边形关系。
// 所有多边形约定按逆时针(CCW)给出边时，内部位于每条有向边的左侧，
// 因此“点到边的有符号距离”为正表示在内部、为负表示在外部、为 0 表示恰在边界上。
//
// 数值策略：跨量级的有限坐标（如宽 1e308、高 1e-100 的矩形）下，
// 任何“先按最大分量归一化再相乘”的做法都会让微小分量下溢为 0、
// 或让乘积溢出为 Infinity，从而把非退化几何误判为退化。
// 因此方向判定、行列式、面积与有符号距离一律改用 BigInt 精确整数运算：
// 每个有限 double 都可精确分解为 整数尾数 × 2^指数，行列式/面积在该表示下
// 精确求值，最后一次性舍入回 double（幅值超界时饱和为 ±Infinity/0，绝不产生 NaN）。

export function sub(a, b) {
  return { x: a.x - b.x, y: a.y - b.y };
}

export function cross(a, b) {
  return a.x * b.y - a.y * b.x;
}

export function pointDistance(a, b) {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  // 先除以较大分量再做 hypot，避免两个超大分量平方和溢出
  const m = Math.max(Math.abs(dx), Math.abs(dy));
  if (m === 0) return 0;
  return Math.hypot(dx / m, dy / m) * m;
}

// ---- 精确算术基础：有限 double ⇔ 整数尾数 × 2^指数 ----

const f64View = new DataView(new ArrayBuffer(8));

/**
 * 把有限 double x 精确分解为 [m, e]：x = m × 2^e（m 为 BigInt）。
 * x 为 0 时返回 [0n, 0]；次正规数同样精确。
 */
function decompose(x) {
  if (x === 0) return [0n, 0];
  f64View.setFloat64(0, x);
  const hi = f64View.getUint32(0);
  const lo = f64View.getUint32(4);
  let mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  const expBits = (hi >>> 20) & 0x7ff;
  let e;
  if (expBits === 0) {
    e = -1074; // 次正规数：无隐含最高位
  } else {
    mant |= 1n << 52n; // 正规数：补上隐含的最高位
    e = expBits - 1075; // 阶码偏置 1023，尾数占 52 位
  }
  return [(hi >>> 31 ? -1n : 1n) * mant, e];
}

/**
 * 两个有限 double 之差的精确值 [m, e]：x - y = m × 2^e。
 * 差超出 double 范围（如 1.8e308 - (-1.8e308)）时依然精确。
 */
function exactDiff(x, y) {
  const [mx, ex] = decompose(x);
  const [my, ey] = decompose(y);
  const e = Math.min(ex, ey);
  return [(mx << BigInt(ex - e)) - (my << BigInt(ey - e)), e];
}

/**
 * 把非零精确值 d × 2^e 规范化为 [h, eh]：h ∈ [2^52, 2^53)，
 * 真实值 ≈ h × 2^eh（截断到 53 位，相对误差 < 2^-52）。d 不得为 0n。
 */
function normalizeExact(d, e) {
  const a = d < 0n ? -d : d;
  const bitLen = a.toString(2).length;
  if (bitLen > 53) return [Number(a >> BigInt(bitLen - 53)), e + bitLen - 53];
  return [Number(a), e];
}

/**
 * 计算 r × 2^e：分段乘 2 的幂，避免 2^e 自身溢出/下溢；
 * 结果超出 double 范围时自然饱和为 Infinity/0。
 */
function scaleByPow2(r, e) {
  while (e > 0 && Number.isFinite(r) && r !== 0) {
    const k = Math.min(e, 1023);
    r *= 2 ** k;
    e -= k;
  }
  while (e < 0 && r !== 0) {
    const k = Math.max(e, -1074);
    r *= 2 ** k;
    e -= k;
  }
  return r;
}

/**
 * 精确值 d × 2^e 舍入为 double：符号始终与 d 一致，
 * 幅值超界时饱和为 ±Infinity，过小则逐渐下溢到 0。
 */
function exactToDouble(d, e) {
  if (d === 0n) return 0;
  const [h, eh] = normalizeExact(d, e);
  const r = scaleByPow2(h, eh);
  return d < 0n ? -r : r;
}

/**
 * 数值稳定的叉积符号：返回 (b-a)×(c-a) 的符号（1/0/-1）。
 * 坐标差与行列式全部用 BigInt 精确求值：即使坐标横跨整个有限 double 范围
 * （如 1e308 与 1e-100 共存、分量相差 1e408 倍），也不会因归一化下溢
 * 或乘积溢出而误判符号。
 */
export function crossSign(a, b, c) {
  const [u1, e1] = exactDiff(b.x, a.x);
  const [v1, f1] = exactDiff(b.y, a.y);
  const [u2, e2] = exactDiff(c.x, a.x);
  const [v2, f2] = exactDiff(c.y, a.y);
  // (b-a)×(c-a) = u1·v2 − v1·u2，两个乘积对齐到共同指数后精确相减
  const g1 = e1 + f2;
  const g2 = f1 + e2;
  const g = Math.min(g1, g2);
  const d = ((u1 * v2) << BigInt(g1 - g)) - ((v1 * u2) << BigInt(g2 - g));
  return d > 0n ? 1 : d < 0n ? -1 : 0;
}

/**
 * 点 p 到过 a、b 的直线的有符号距离（沿 a->b 方向，左侧为正）。
 * 分子（行列式）与分母（边长）都按 尾数×2^指数 的精确形式求值后再相除，
 * 保证坐标在整个有限 double 范围内、任意长宽比下都返回正确量级的有限值，
 * 而不会因中间量下溢（如 5e-101/1e308 → 0）把正裕量错误归零。
 */
export function lineSide(p, a, b) {
  if (a.x === b.x && a.y === b.y) {
    // a、b 重合：退化为到点 a 的距离
    const qx = p.x - a.x;
    const qy = p.y - a.y;
    const m = Math.max(Math.abs(qx), Math.abs(qy));
    return m === 0 ? 0 : Math.hypot(qx / m, qy / m) * m;
  }
  const [dx, edx] = exactDiff(b.x, a.x);
  const [dy, edy] = exactDiff(b.y, a.y);
  const [qx, eqx] = exactDiff(p.x, a.x);
  const [qy, eqy] = exactDiff(p.y, a.y);
  // 有符号距离 = (dx·qy − dy·qx) / |d|，分子精确
  const g1 = edx + eqy;
  const g2 = edy + eqx;
  const g = Math.min(g1, g2);
  const det = ((dx * qy) << BigInt(g1 - g)) - ((dy * qx) << BigInt(g2 - g));
  if (det === 0n) return 0;
  // |d|² = dx² + dy² 同样精确，再对精确值开方，避免 hypot 的中间溢出
  const h1 = 2 * edx;
  const h2 = 2 * edy;
  const h = Math.min(h1, h2);
  const sumSq = ((dx * dx) << BigInt(h1 - h)) + ((dy * dy) << BigInt(h2 - h));
  let [h0, eh] = normalizeExact(sumSq, h); // |d|² ≈ h0 × 2^eh
  if (eh % 2 !== 0) {
    h0 *= 2;
    eh -= 1;
  }
  const lenMant = Math.sqrt(h0);
  const lenExp = eh / 2;
  // 商 = det×2^g / (lenMant×2^lenExp)：尾数相除（量级相近、不会溢出）后再分段放大
  const [detMant, detExp] = normalizeExact(det, g);
  const r = scaleByPow2(detMant / lenMant, detExp - lenExp);
  return det < 0n ? -r : r;
}

/**
 * 点到线段的（无符号）距离，用于退化凸包（共线/重合）时的证据量化。
 * 投影参数按统一尺度归一化计算，避免超大坐标下相乘溢出。
 */
export function distanceToSegment(p, a, b) {
  const dx = b.x - a.x;
  const dy = b.y - a.y;
  const px = p.x - a.x;
  const py = p.y - a.y;
  if (dx === 0 && dy === 0) {
    const m = Math.max(Math.abs(px), Math.abs(py));
    return m === 0 ? 0 : Math.hypot(px / m, py / m) * m;
  }
  const s = Math.max(Math.abs(dx), Math.abs(dy), Math.abs(px), Math.abs(py));
  const dxs = dx / s;
  const dys = dy / s;
  const len2 = dxs * dxs + dys * dys;
  if (len2 === 0) {
    return Math.hypot(px / s, py / s) * s;
  }
  let t = ((px / s) * dxs + (py / s) * dys) / len2;
  t = Math.max(0, Math.min(1, t));
  const rx = px / s - t * dxs;
  const ry = py / s - t * dys;
  return Math.hypot(rx, ry) * s;
}

export function polygonSignedArea(poly) {
  // shoelace 和 Σ(a.x·b.y − b.x·a.y) 用 BigInt 精确累加：
  // 若先按最大坐标归一化，1e-100 量级的坐标相对 1e308 会下溢为 0，
  // 会把面积为 1e208 的窄条矩形误判成退化（面积为 0）。
  let minExp = Infinity;
  const terms = [];
  const push = (x, y, sign) => {
    const [mx, ex] = decompose(x);
    const [my, ey] = decompose(y);
    if (mx === 0n || my === 0n) return;
    const e = ex + ey;
    if (e < minExp) minExp = e;
    terms.push([sign * mx * my, e]);
  };
  for (let i = 0; i < poly.length; i++) {
    const a = poly[i];
    const b = poly[(i + 1) % poly.length];
    push(a.x, b.y, 1n);
    push(b.x, a.y, -1n);
  }
  if (terms.length === 0) return 0;
  let sum = 0n;
  for (const [m, e] of terms) sum += m << BigInt(e - minExp);
  // 面积 = 总和 / 2 → 指数减 1；真实面积超过 double 上限时饱和为 ±Infinity，
  // 仍保留正确符号；退化（共线）情形 sum 为 0，返回 0。
  return exactToDouble(sum, minExp - 1);
}

/**
 * 保证多边形为逆时针方向。
 */
export function ensureCCW(poly) {
  return polygonSignedArea(poly) < 0 ? poly.slice().reverse() : poly.slice();
}

/**
 * 点到（逆时针）多边形边界的有符号距离：
 * 内部为正（到最近边的距离），外部为负。
 */
export function polygonSignedDistance(polyCCW, p) {
  let minAbs = Infinity;
  for (let i = 0; i < polyCCW.length; i++) {
    minAbs = Math.min(minAbs, lineSide(p, polyCCW[i], polyCCW[(i + 1) % polyCCW.length]));
  }
  return minAbs;
}

/**
 * 射线法点是否在多边形内（含边界）。
 * 边界共线判定与射线交点均用不产生超大乘积的形式，
 * 保证在整个有限 double 坐标范围内结果可靠。
 */
export function pointInPolygon(poly, p) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const a = poly[i];
    const b = poly[j];
    const onBoundary =
      crossSign(a, b, p) === 0 &&
      Math.min(a.x, b.x) <= p.x && p.x <= Math.max(a.x, b.x) &&
      Math.min(a.y, b.y) <= p.y && p.y <= Math.max(a.y, b.y);
    if (onBoundary) return true;
    if (a.y > p.y !== b.y > p.y) {
      // 交点参数 t∈(0,1)，t*(b.x-a.x) 的量级不超过边长，避免两个超大坐标相乘溢出
      const t = (p.y - a.y) / (b.y - a.y);
      const xIntersect = a.x + t * (b.x - a.x);
      if (p.x < xIntersect) inside = !inside;
    }
  }
  return inside;
}

function dedupe(points) {
  const seen = new Set();
  const out = [];
  for (const p of points) {
    const key = `${p.x},${p.y}`;
    if (!seen.has(key)) {
      seen.add(key);
      out.push(p);
    }
  }
  return out;
}

/**
 * Andrew 单调链凸包，返回逆时针顶点序列，剔除共线中间点。
 * 点数不足时原样返回（去重后）。
 */
export function convexHull(points) {
  const pts = dedupe(points).sort((a, b) => a.x - b.x || a.y - b.y);
  if (pts.length <= 2) return pts;

  const lower = [];
  for (const p of pts) {
    while (lower.length >= 2 && crossSign(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
      lower.pop();
    }
    lower.push(p);
  }
  const upper = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    const p = pts[i];
    while (upper.length >= 2 && crossSign(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
      upper.pop();
    }
    upper.push(p);
  }
  lower.pop();
  upper.pop();
  return lower.concat(upper);
}

/**
 * 点到逆时针凸包边界的最小有符号距离（各边有符号距离的最小值）。
 * 凸包退化（不足 3 个顶点）时返回 -Infinity。
 */
export function convexMargin(hullCCW, p) {
  if (hullCCW.length < 3) return -Infinity;
  let m = Infinity;
  for (let i = 0; i < hullCCW.length; i++) {
    m = Math.min(m, lineSide(p, hullCCW[i], hullCCW[(i + 1) % hullCCW.length]));
  }
  return m;
}
