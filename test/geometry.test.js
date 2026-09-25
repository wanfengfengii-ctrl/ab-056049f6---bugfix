import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  convexHull,
  ensureCCW,
  convexMargin,
  lineSide,
  pointInPolygon,
  polygonSignedArea,
  polygonSignedDistance,
  pointDistance,
  distanceToSegment,
  crossSign,
} from '../server/geometry.js';

test('凸包按逆时针返回并剔除共线点', () => {
  const hull = ensureCCW(convexHull([
    { x: 0, y: 0 }, { x: 2, y: 0 }, { x: 2, y: 2 },
    { x: 0, y: 2 }, { x: 1, y: 1 },
  ]));
  assert.equal(hull.length, 4);
  assert.ok(polygonSignedArea(hull) > 0);
});

test('内部点有符号距离为正，外部为负', () => {
  const hull = ensureCCW([
    { x: -5, y: -5 }, { x: 5, y: -5 }, { x: 5, y: 5 }, { x: -5, y: 5 },
  ]);
  assert.ok(convexMargin(hull, { x: 0, y: 0 }) > 0);
  assert.equal(convexMargin(hull, { x: 0, y: 0 }).toFixed(3), '5.000');
  assert.ok(convexMargin(hull, { x: 0, y: 6 }) < 0);
});

test('恰在边上的点有符号距离为 0（不满足严格在内）', () => {
  const hull = ensureCCW([
    { x: -5, y: -5 }, { x: 5, y: -5 }, { x: 5, y: 5 }, { x: -5, y: 5 },
  ]);
  assert.equal(convexMargin(hull, { x: 0, y: 5 }), 0);
});

test('lineSide 左侧为正', () => {
  assert.ok(lineSide({ x: 0, y: 1 }, { x: -1, y: 0 }, { x: 1, y: 0 }) > 0);
  assert.ok(lineSide({ x: 0, y: -1 }, { x: -1, y: 0 }, { x: 1, y: 0 }) < 0);
});

test('pointInPolygon 含边界', () => {
  const poly = [{ x: 0, y: 0 }, { x: 4, y: 0 }, { x: 4, y: 4 }, { x: 0, y: 4 }];
  assert.equal(pointInPolygon(poly, { x: 2, y: 2 }), true);
  assert.equal(pointInPolygon(poly, { x: 4, y: 2 }), true);
  assert.equal(pointInPolygon(poly, { x: 5, y: 2 }), false);
});

// ---- 超大有限坐标（1e307 量级接近 double 上界 1.8e308）下的数值稳定性 ----

const H = 1e307;
const bigSquareCCW = [
  { x: -H, y: -H }, { x: H, y: -H }, { x: H, y: H }, { x: -H, y: H },
];

test('超大坐标：lineSide 不溢出，原点到各边有符号距离约为 1e307', () => {
  assert.ok(Number.isFinite(lineSide({ x: 0, y: 0 }, bigSquareCCW[0], bigSquareCCW[1])));
  assert.ok(Math.abs(lineSide({ x: 0, y: 0 }, bigSquareCCW[0], bigSquareCCW[1]) - H) / H < 1e-12);
  // 顶边走向为右→左，其左侧（内部）在下方；外侧点 y=H+1e291 的有符号距离应为负且有限
  const outside = lineSide({ x: 0, y: H + 1e291 }, bigSquareCCW[2], bigSquareCCW[3]);
  assert.ok(Number.isFinite(outside) && outside < 0, `外侧距离异常：${outside}`);
});

test('超大坐标：convexMargin 返回有限值而非 null/Infinity/NaN', () => {
  const m = convexMargin(bigSquareCCW, { x: 0, y: 0 });
  assert.ok(Number.isFinite(m));
  assert.ok(Math.abs(m - H) / H < 1e-12);
});

test('超大坐标：凸包正确识别 4 个角且面积符号为正（不因溢出变成 NaN）', () => {
  const hull = convexHull(bigSquareCCW);
  assert.equal(hull.length, 4);
  // 面积真值约 4e614 超出 double，归一化后返回 +Infinity，但绝不能是 NaN/负值
  const area = polygonSignedArea(ensureCCW(bigSquareCCW));
  assert.ok(!Number.isNaN(area) && area > 0, `面积符号异常：${area}`);
});

test('超大坐标：原点在超大方形内，外部点判定为不在', () => {
  assert.equal(pointInPolygon(bigSquareCCW, { x: 0, y: 0 }), true);
  assert.equal(pointInPolygon(bigSquareCCW, { x: 1.5e307, y: 0 }), false);
});

test('超大坐标：pointDistance 不溢出', () => {
  const d = pointDistance({ x: -H, y: -H }, { x: H, y: H });
  assert.ok(Number.isFinite(d));
  assert.ok(Math.abs(d - 2 * H * Math.SQRT2) / (2 * H * Math.SQRT2) < 1e-12);
});

test('超大坐标：distanceToSegment 内/外投影均有限', () => {
  const a = { x: -H, y: -H };
  const b = { x: H, y: -H };
  const inside = distanceToSegment({ x: 0, y: 0 }, a, b);
  const beyond = distanceToSegment({ x: 0, y: -H - 5e306 }, a, b);
  assert.ok(Number.isFinite(inside) && inside > 0);
  assert.ok(Number.isFinite(beyond) && beyond > 0);
  assert.equal(Number.isNaN(distanceToSegment(a, a, b)), false);
});

test('超大坐标：crossSign 对明显左/右/共线给出正确符号', () => {
  assert.equal(crossSign({ x: -H, y: 0 }, { x: H, y: 0 }, { x: 0, y: 1 }), 1);
  assert.equal(crossSign({ x: -H, y: 0 }, { x: H, y: 0 }, { x: 0, y: -1 }), -1);
  assert.equal(crossSign({ x: -H, y: 0 }, { x: H, y: 0 }, { x: 0, y: 0 }), 0);
});

test('超大坐标：多边形有符号距离返回有限值', () => {
  const d = polygonSignedDistance(bigSquareCCW, { x: 0, y: 0 });
  assert.ok(Number.isFinite(d));
  assert.ok(d > 0);
});

// ---- 极端有限长宽比（宽 1e308、高 1e-100）下不得误判为退化 ----

const EW = 1e308;
const EH = 1e-100;
const extremeRect = [
  { x: 0, y: 0 }, { x: EW, y: 0 }, { x: EW, y: EH }, { x: 0, y: EH },
];

test('极端长宽比：crossSign 识别四个角点非共线（真值叉积 1e208）', () => {
  // 前三点：(0,0)->(1e308,0)->(1e308,1e-100)，叉积真值为正且巨大
  assert.equal(crossSign(extremeRect[0], extremeRect[1], extremeRect[2]), 1);
  // 反向拐角
  assert.equal(crossSign(extremeRect[1], extremeRect[2], extremeRect[3]), 1);
  // 真正共线（细长方向）判 0
  assert.equal(
    crossSign({ x: 0, y: 0 }, { x: EW, y: 0 }, { x: EW / 2, y: 0 }),
    0
  );
});

test('极端长宽比：convexHull 返回全部 4 个顶点', () => {
  const hull = convexHull(extremeRect);
  assert.equal(hull.length, 4, `hull=${JSON.stringify(hull)}`);
  assert.ok(polygonSignedArea(ensureCCW(hull)) > 0);
});

test('极端长宽比：polygonSignedArea 约为 1e208（有限正数，非 0/Infinity/NaN）', () => {
  const area = polygonSignedArea(extremeRect);
  assert.ok(Number.isFinite(area), `area 必须有限，实际 ${area}`);
  assert.ok(Math.abs(area - 1e208) / 1e208 < 1e-12, `area=${area}`);
});

test('极端长宽比：内部点到水平边的有符号距离约为 5e-101，到竖直边约 5e307', () => {
  const ccw = ensureCCW(extremeRect);
  const p = { x: 5e307, y: 5e-101 };
  // convexMargin 取各边有符号距离最小值，由上下水平边锁定为 5e-101
  const m = convexMargin(ccw, p);
  assert.ok(Number.isFinite(m), `margin 必须有限，实际 ${m}`);
  assert.ok(m > 0, '内部点裕量必须严格为正');
  assert.ok(Math.abs(m - 5e-101) / 5e-101 < 1e-12, `margin=${m}`);
});

test('极端长宽比：pointInPolygon 对矩形内部/边上/外部判定正确', () => {
  assert.equal(pointInPolygon(extremeRect, { x: 5e307, y: 5e-101 }), true);
  // 落在水平边上 => 含边界 => true
  assert.equal(pointInPolygon(extremeRect, { x: 5e307, y: 0 }), true);
  // y 方向越界（与 1e-100 相比可分辨）
  assert.equal(pointInPolygon(extremeRect, { x: 5e307, y: 2e-100 }), false);
});
