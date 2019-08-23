// modules are defined as an array
// [ module function, map of requires ]
//
// map of requires is short require name -> numeric require
//
// anything defined in a previous bundle is accessed via the
// orig method which is the require for previous bundles
parcelRequire = (function (modules, cache, entry, globalName) {
  // Save the require from previous bundle to this closure if any
  var previousRequire = typeof parcelRequire === 'function' && parcelRequire;
  var nodeRequire = typeof require === 'function' && require;

  function newRequire(name, jumped) {
    if (!cache[name]) {
      if (!modules[name]) {
        // if we cannot find the module within our internal map or
        // cache jump to the current global require ie. the last bundle
        // that was added to the page.
        var currentRequire = typeof parcelRequire === 'function' && parcelRequire;
        if (!jumped && currentRequire) {
          return currentRequire(name, true);
        }

        // If there are other bundles on this page the require from the
        // previous one is saved to 'previousRequire'. Repeat this as
        // many times as there are bundles until the module is found or
        // we exhaust the require chain.
        if (previousRequire) {
          return previousRequire(name, true);
        }

        // Try the node require function if it exists.
        if (nodeRequire && typeof name === 'string') {
          return nodeRequire(name);
        }

        var err = new Error('Cannot find module \'' + name + '\'');
        err.code = 'MODULE_NOT_FOUND';
        throw err;
      }

      localRequire.resolve = resolve;
      localRequire.cache = {};

      var module = cache[name] = new newRequire.Module(name);

      modules[name][0].call(module.exports, localRequire, module, module.exports, this);
    }

    return cache[name].exports;

    function localRequire(x){
      return newRequire(localRequire.resolve(x));
    }

    function resolve(x){
      return modules[name][1][x] || x;
    }
  }

  function Module(moduleName) {
    this.id = moduleName;
    this.bundle = newRequire;
    this.exports = {};
  }

  newRequire.isParcelRequire = true;
  newRequire.Module = Module;
  newRequire.modules = modules;
  newRequire.cache = cache;
  newRequire.parent = previousRequire;
  newRequire.register = function (id, exports) {
    modules[id] = [function (require, module) {
      module.exports = exports;
    }, {}];
  };

  var error;
  for (var i = 0; i < entry.length; i++) {
    try {
      newRequire(entry[i]);
    } catch (e) {
      // Save first error but execute all entries
      if (!error) {
        error = e;
      }
    }
  }

  if (entry.length) {
    // Expose entry point to Node, AMD or browser globals
    // Based on https://github.com/ForbesLindesay/umd/blob/master/template.js
    var mainExports = newRequire(entry[entry.length - 1]);

    // CommonJS
    if (typeof exports === "object" && typeof module !== "undefined") {
      module.exports = mainExports;

    // RequireJS
    } else if (typeof define === "function" && define.amd) {
     define(function () {
       return mainExports;
     });

    // <script>
    } else if (globalName) {
      this[globalName] = mainExports;
    }
  }

  // Override the current require with this new one
  parcelRequire = newRequire;

  if (error) {
    // throw error from earlier, _after updating parcelRequire_
    throw error;
  }

  return newRequire;
})({"math/vec2.ts":[function(require,module,exports) {
"use strict";

exports.__esModule = true;

var Vec2 =
/** @class */
function () {
  function Vec2(x, y) {
    this.x = x;
    this.y = y;
  }

  return Vec2;
}();

exports.Vec2 = Vec2;
},{}],"../node_modules/delaunator/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;
const EPSILON = Math.pow(2, -52);
const EDGE_STACK = new Uint32Array(512);

class Delaunator {
  static from(points, getX = defaultGetX, getY = defaultGetY) {
    const n = points.length;
    const coords = new Float64Array(n * 2);

    for (let i = 0; i < n; i++) {
      const p = points[i];
      coords[2 * i] = getX(p);
      coords[2 * i + 1] = getY(p);
    }

    return new Delaunator(coords);
  }

  constructor(coords) {
    const n = coords.length >> 1;
    if (n > 0 && typeof coords[0] !== 'number') throw new Error('Expected coords to contain numbers.');
    this.coords = coords; // arrays that will store the triangulation graph

    const maxTriangles = Math.max(2 * n - 5, 0);
    this._triangles = new Uint32Array(maxTriangles * 3);
    this._halfedges = new Int32Array(maxTriangles * 3); // temporary arrays for tracking the edges of the advancing convex hull

    this._hashSize = Math.ceil(Math.sqrt(n));
    this._hullPrev = new Uint32Array(n); // edge to prev edge

    this._hullNext = new Uint32Array(n); // edge to next edge

    this._hullTri = new Uint32Array(n); // edge to adjacent triangle

    this._hullHash = new Int32Array(this._hashSize).fill(-1); // angular edge hash
    // temporary arrays for sorting points

    this._ids = new Uint32Array(n);
    this._dists = new Float64Array(n);
    this.update();
  }

  update() {
    const {
      coords,
      _hullPrev: hullPrev,
      _hullNext: hullNext,
      _hullTri: hullTri,
      _hullHash: hullHash
    } = this;
    const n = coords.length >> 1; // populate an array of point indices; calculate input data bbox

    let minX = Infinity;
    let minY = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;

    for (let i = 0; i < n; i++) {
      const x = coords[2 * i];
      const y = coords[2 * i + 1];
      if (x < minX) minX = x;
      if (y < minY) minY = y;
      if (x > maxX) maxX = x;
      if (y > maxY) maxY = y;
      this._ids[i] = i;
    }

    const cx = (minX + maxX) / 2;
    const cy = (minY + maxY) / 2;
    let minDist = Infinity;
    let i0, i1, i2; // pick a seed point close to the center

    for (let i = 0; i < n; i++) {
      const d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);

      if (d < minDist) {
        i0 = i;
        minDist = d;
      }
    }

    const i0x = coords[2 * i0];
    const i0y = coords[2 * i0 + 1];
    minDist = Infinity; // find the point closest to the seed

    for (let i = 0; i < n; i++) {
      if (i === i0) continue;
      const d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);

      if (d < minDist && d > 0) {
        i1 = i;
        minDist = d;
      }
    }

    let i1x = coords[2 * i1];
    let i1y = coords[2 * i1 + 1];
    let minRadius = Infinity; // find the third point which forms the smallest circumcircle with the first two

    for (let i = 0; i < n; i++) {
      if (i === i0 || i === i1) continue;
      const r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);

      if (r < minRadius) {
        i2 = i;
        minRadius = r;
      }
    }

    let i2x = coords[2 * i2];
    let i2y = coords[2 * i2 + 1];

    if (minRadius === Infinity) {
      // order collinear points by dx (or dy if all x are identical)
      // and return the list as a hull
      for (let i = 0; i < n; i++) {
        this._dists[i] = coords[2 * i] - coords[0] || coords[2 * i + 1] - coords[1];
      }

      quicksort(this._ids, this._dists, 0, n - 1);
      const hull = new Uint32Array(n);
      let j = 0;

      for (let i = 0, d0 = -Infinity; i < n; i++) {
        const id = this._ids[i];

        if (this._dists[id] > d0) {
          hull[j++] = id;
          d0 = this._dists[id];
        }
      }

      this.hull = hull.subarray(0, j);
      this.triangles = new Uint32Array(0);
      this.halfedges = new Uint32Array(0);
      return;
    } // swap the order of the seed points for counter-clockwise orientation


    if (orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
      const i = i1;
      const x = i1x;
      const y = i1y;
      i1 = i2;
      i1x = i2x;
      i1y = i2y;
      i2 = i;
      i2x = x;
      i2y = y;
    }

    const center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
    this._cx = center.x;
    this._cy = center.y;

    for (let i = 0; i < n; i++) {
      this._dists[i] = dist(coords[2 * i], coords[2 * i + 1], center.x, center.y);
    } // sort the points by distance from the seed triangle circumcenter


    quicksort(this._ids, this._dists, 0, n - 1); // set up the seed triangle as the starting hull

    this._hullStart = i0;
    let hullSize = 3;
    hullNext[i0] = hullPrev[i2] = i1;
    hullNext[i1] = hullPrev[i0] = i2;
    hullNext[i2] = hullPrev[i1] = i0;
    hullTri[i0] = 0;
    hullTri[i1] = 1;
    hullTri[i2] = 2;
    hullHash.fill(-1);
    hullHash[this._hashKey(i0x, i0y)] = i0;
    hullHash[this._hashKey(i1x, i1y)] = i1;
    hullHash[this._hashKey(i2x, i2y)] = i2;
    this.trianglesLen = 0;

    this._addTriangle(i0, i1, i2, -1, -1, -1);

    for (let k = 0, xp, yp; k < this._ids.length; k++) {
      const i = this._ids[k];
      const x = coords[2 * i];
      const y = coords[2 * i + 1]; // skip near-duplicate points

      if (k > 0 && Math.abs(x - xp) <= EPSILON && Math.abs(y - yp) <= EPSILON) continue;
      xp = x;
      yp = y; // skip seed triangle points

      if (i === i0 || i === i1 || i === i2) continue; // find a visible edge on the convex hull using edge hash

      let start = 0;

      for (let j = 0, key = this._hashKey(x, y); j < this._hashSize; j++) {
        start = hullHash[(key + j) % this._hashSize];
        if (start !== -1 && start !== hullNext[start]) break;
      }

      start = hullPrev[start];
      let e = start,
          q;

      while (q = hullNext[e], !orient(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1])) {
        e = q;

        if (e === start) {
          e = -1;
          break;
        }
      }

      if (e === -1) continue; // likely a near-duplicate point; skip it
      // add the first triangle from the point

      let t = this._addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]); // recursively flip triangles from the point until they satisfy the Delaunay condition


      hullTri[i] = this._legalize(t + 2);
      hullTri[e] = t; // keep track of boundary triangles on the hull

      hullSize++; // walk forward through the hull, adding more triangles and flipping recursively

      let n = hullNext[e];

      while (q = hullNext[n], orient(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1])) {
        t = this._addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
        hullTri[i] = this._legalize(t + 2);
        hullNext[n] = n; // mark as removed

        hullSize--;
        n = q;
      } // walk backward from the other side, adding more triangles and flipping


      if (e === start) {
        while (q = hullPrev[e], orient(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1])) {
          t = this._addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);

          this._legalize(t + 2);

          hullTri[q] = t;
          hullNext[e] = e; // mark as removed

          hullSize--;
          e = q;
        }
      } // update the hull indices


      this._hullStart = hullPrev[i] = e;
      hullNext[e] = hullPrev[n] = i;
      hullNext[i] = n; // save the two new edges in the hash table

      hullHash[this._hashKey(x, y)] = i;
      hullHash[this._hashKey(coords[2 * e], coords[2 * e + 1])] = e;
    }

    this.hull = new Uint32Array(hullSize);

    for (let i = 0, e = this._hullStart; i < hullSize; i++) {
      this.hull[i] = e;
      e = hullNext[e];
    } // trim typed triangle mesh arrays


    this.triangles = this._triangles.subarray(0, this.trianglesLen);
    this.halfedges = this._halfedges.subarray(0, this.trianglesLen);
  }

  _hashKey(x, y) {
    return Math.floor(pseudoAngle(x - this._cx, y - this._cy) * this._hashSize) % this._hashSize;
  }

  _legalize(a) {
    const {
      _triangles: triangles,
      _halfedges: halfedges,
      coords
    } = this;
    let i = 0;
    let ar = 0; // recursion eliminated with a fixed-size stack

    while (true) {
      const b = halfedges[a];
      /* if the pair of triangles doesn't satisfy the Delaunay condition
       * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
       * then do the same check/flip recursively for the new pair of triangles
       *
       *           pl                    pl
       *          /||\                  /  \
       *       al/ || \bl            al/    \a
       *        /  ||  \              /      \
       *       /  a||b  \    flip    /___ar___\
       *     p0\   ||   /p1   =>   p0\---bl---/p1
       *        \  ||  /              \      /
       *       ar\ || /br             b\    /br
       *          \||/                  \  /
       *           pr                    pr
       */

      const a0 = a - a % 3;
      ar = a0 + (a + 2) % 3;

      if (b === -1) {
        // convex hull edge
        if (i === 0) break;
        a = EDGE_STACK[--i];
        continue;
      }

      const b0 = b - b % 3;
      const al = a0 + (a + 1) % 3;
      const bl = b0 + (b + 2) % 3;
      const p0 = triangles[ar];
      const pr = triangles[a];
      const pl = triangles[al];
      const p1 = triangles[bl];
      const illegal = inCircle(coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl], coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]);

      if (illegal) {
        triangles[a] = p1;
        triangles[b] = p0;
        const hbl = halfedges[bl]; // edge swapped on the other side of the hull (rare); fix the halfedge reference

        if (hbl === -1) {
          let e = this._hullStart;

          do {
            if (this._hullTri[e] === bl) {
              this._hullTri[e] = a;
              break;
            }

            e = this._hullPrev[e];
          } while (e !== this._hullStart);
        }

        this._link(a, hbl);

        this._link(b, halfedges[ar]);

        this._link(ar, bl);

        const br = b0 + (b + 1) % 3; // don't worry about hitting the cap: it can only happen on extremely degenerate input

        if (i < EDGE_STACK.length) {
          EDGE_STACK[i++] = br;
        }
      } else {
        if (i === 0) break;
        a = EDGE_STACK[--i];
      }
    }

    return ar;
  }

  _link(a, b) {
    this._halfedges[a] = b;
    if (b !== -1) this._halfedges[b] = a;
  } // add a new triangle given vertex indices and adjacent half-edge ids


  _addTriangle(i0, i1, i2, a, b, c) {
    const t = this.trianglesLen;
    this._triangles[t] = i0;
    this._triangles[t + 1] = i1;
    this._triangles[t + 2] = i2;

    this._link(t, a);

    this._link(t + 1, b);

    this._link(t + 2, c);

    this.trianglesLen += 3;
    return t;
  }

} // monotonically increases with real angle, but doesn't need expensive trigonometry


exports.default = Delaunator;

function pseudoAngle(dx, dy) {
  const p = dx / (Math.abs(dx) + Math.abs(dy));
  return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
}

function dist(ax, ay, bx, by) {
  const dx = ax - bx;
  const dy = ay - by;
  return dx * dx + dy * dy;
}

function orient(px, py, qx, qy, rx, ry) {
  return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0;
}

function inCircle(ax, ay, bx, by, cx, cy, px, py) {
  const dx = ax - px;
  const dy = ay - py;
  const ex = bx - px;
  const ey = by - py;
  const fx = cx - px;
  const fy = cy - py;
  const ap = dx * dx + dy * dy;
  const bp = ex * ex + ey * ey;
  const cp = fx * fx + fy * fy;
  return dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0;
}

function circumradius(ax, ay, bx, by, cx, cy) {
  const dx = bx - ax;
  const dy = by - ay;
  const ex = cx - ax;
  const ey = cy - ay;
  const bl = dx * dx + dy * dy;
  const cl = ex * ex + ey * ey;
  const d = 0.5 / (dx * ey - dy * ex);
  const x = (ey * bl - dy * cl) * d;
  const y = (dx * cl - ex * bl) * d;
  return x * x + y * y;
}

function circumcenter(ax, ay, bx, by, cx, cy) {
  const dx = bx - ax;
  const dy = by - ay;
  const ex = cx - ax;
  const ey = cy - ay;
  const bl = dx * dx + dy * dy;
  const cl = ex * ex + ey * ey;
  const d = 0.5 / (dx * ey - dy * ex);
  const x = ax + (ey * bl - dy * cl) * d;
  const y = ay + (dx * cl - ex * bl) * d;
  return {
    x,
    y
  };
}

function quicksort(ids, dists, left, right) {
  if (right - left <= 20) {
    for (let i = left + 1; i <= right; i++) {
      const temp = ids[i];
      const tempDist = dists[temp];
      let j = i - 1;

      while (j >= left && dists[ids[j]] > tempDist) ids[j + 1] = ids[j--];

      ids[j + 1] = temp;
    }
  } else {
    const median = left + right >> 1;
    let i = left + 1;
    let j = right;
    swap(ids, median, i);
    if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
    if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
    if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);
    const temp = ids[i];
    const tempDist = dists[temp];

    while (true) {
      do i++; while (dists[ids[i]] < tempDist);

      do j--; while (dists[ids[j]] > tempDist);

      if (j < i) break;
      swap(ids, i, j);
    }

    ids[left + 1] = ids[j];
    ids[j] = temp;

    if (right - i + 1 >= j - left) {
      quicksort(ids, dists, i, right);
      quicksort(ids, dists, left, j - 1);
    } else {
      quicksort(ids, dists, left, j - 1);
      quicksort(ids, dists, i, right);
    }
  }
}

function swap(arr, i, j) {
  const tmp = arr[i];
  arr[i] = arr[j];
  arr[j] = tmp;
}

function defaultGetX(p) {
  return p[0];
}

function defaultGetY(p) {
  return p[1];
}
},{}],"../node_modules/d3-delaunay/src/path.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;
const epsilon = 1e-6;

class Path {
  constructor() {
    this._x0 = this._y0 = // start of current subpath
    this._x1 = this._y1 = null; // end of current subpath

    this._ = "";
  }

  moveTo(x, y) {
    this._ += `M${this._x0 = this._x1 = +x},${this._y0 = this._y1 = +y}`;
  }

  closePath() {
    if (this._x1 !== null) {
      this._x1 = this._x0, this._y1 = this._y0;
      this._ += "Z";
    }
  }

  lineTo(x, y) {
    this._ += `L${this._x1 = +x},${this._y1 = +y}`;
  }

  arc(x, y, r) {
    x = +x, y = +y, r = +r;
    const x0 = x + r;
    const y0 = y;
    if (r < 0) throw new Error("negative radius");
    if (this._x1 === null) this._ += `M${x0},${y0}`;else if (Math.abs(this._x1 - x0) > epsilon || Math.abs(this._y1 - y0) > epsilon) this._ += "L" + x0 + "," + y0;
    if (!r) return;
    this._ += `A${r},${r},0,1,1,${x - r},${y}A${r},${r},0,1,1,${this._x1 = x0},${this._y1 = y0}`;
  }

  rect(x, y, w, h) {
    this._ += `M${this._x0 = this._x1 = +x},${this._y0 = this._y1 = +y}h${+w}v${+h}h${-w}Z`;
  }

  value() {
    return this._ || null;
  }

}

exports.default = Path;
},{}],"../node_modules/d3-delaunay/src/polygon.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

class Polygon {
  constructor() {
    this._ = [];
  }

  moveTo(x, y) {
    this._.push([x, y]);
  }

  closePath() {
    this._.push(this._[0].slice());
  }

  lineTo(x, y) {
    this._.push([x, y]);
  }

  value() {
    return this._.length ? this._ : null;
  }

}

exports.default = Polygon;
},{}],"../node_modules/d3-delaunay/src/voronoi.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _path = _interopRequireDefault(require("./path.js"));

var _polygon = _interopRequireDefault(require("./polygon.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

class Voronoi {
  constructor(delaunay, [xmin, ymin, xmax, ymax] = [0, 0, 960, 500]) {
    if (!((xmax = +xmax) >= (xmin = +xmin)) || !((ymax = +ymax) >= (ymin = +ymin))) throw new Error("invalid bounds");
    this.delaunay = delaunay;
    this._circumcenters = new Float64Array(delaunay.points.length * 2);
    this.vectors = new Float64Array(delaunay.points.length * 2);
    this.xmax = xmax, this.xmin = xmin;
    this.ymax = ymax, this.ymin = ymin;

    this._init();
  }

  update() {
    this.delaunay.update();

    this._init();

    return this;
  }

  _init() {
    const {
      delaunay: {
        points,
        hull,
        triangles
      },
      vectors
    } = this; // Compute circumcenters.

    const circumcenters = this.circumcenters = this._circumcenters.subarray(0, triangles.length / 3 * 2);

    for (let i = 0, j = 0, n = triangles.length; i < n; i += 3, j += 2) {
      const t1 = triangles[i] * 2;
      const t2 = triangles[i + 1] * 2;
      const t3 = triangles[i + 2] * 2;
      const x1 = points[t1];
      const y1 = points[t1 + 1];
      const x2 = points[t2];
      const y2 = points[t2 + 1];
      const x3 = points[t3];
      const y3 = points[t3 + 1];
      const a2 = x1 - x2;
      const a3 = x1 - x3;
      const b2 = y1 - y2;
      const b3 = y1 - y3;
      const d1 = x1 * x1 + y1 * y1;
      const d2 = d1 - x2 * x2 - y2 * y2;
      const d3 = d1 - x3 * x3 - y3 * y3;
      const ab = (a3 * b2 - a2 * b3) * 2; // degenerate case (2 points)

      if (!ab) {
        circumcenters[j] = (x1 + x3) / 2 + 1e8 * b3;
        circumcenters[j + 1] = (y1 + y3) / 2 - 1e8 * a3;
      } else {
        circumcenters[j] = (b2 * d3 - b3 * d2) / ab;
        circumcenters[j + 1] = (a3 * d2 - a2 * d3) / ab;
      }
    } // Compute exterior cell rays.


    let h = hull[hull.length - 1];
    let p0,
        p1 = h * 4;
    let x0,
        x1 = points[2 * h];
    let y0,
        y1 = points[2 * h + 1];
    vectors.fill(0);

    for (let i = 0; i < hull.length; ++i) {
      h = hull[i];
      p0 = p1, x0 = x1, y0 = y1;
      p1 = h * 4, x1 = points[2 * h], y1 = points[2 * h + 1];
      vectors[p0 + 2] = vectors[p1] = y0 - y1;
      vectors[p0 + 3] = vectors[p1 + 1] = x1 - x0;
    }
  }

  render(context) {
    const buffer = context == null ? context = new _path.default() : undefined;
    const {
      delaunay: {
        halfedges,
        inedges,
        hull
      },
      circumcenters,
      vectors
    } = this;
    if (hull.length <= 1) return null;

    for (let i = 0, n = halfedges.length; i < n; ++i) {
      const j = halfedges[i];
      if (j < i) continue;
      const ti = Math.floor(i / 3) * 2;
      const tj = Math.floor(j / 3) * 2;
      const xi = circumcenters[ti];
      const yi = circumcenters[ti + 1];
      const xj = circumcenters[tj];
      const yj = circumcenters[tj + 1];

      this._renderSegment(xi, yi, xj, yj, context);
    }

    let h0,
        h1 = hull[hull.length - 1];

    for (let i = 0; i < hull.length; ++i) {
      h0 = h1, h1 = hull[i];
      const t = Math.floor(inedges[h1] / 3) * 2;
      const x = circumcenters[t];
      const y = circumcenters[t + 1];
      const v = h0 * 4;

      const p = this._project(x, y, vectors[v + 2], vectors[v + 3]);

      if (p) this._renderSegment(x, y, p[0], p[1], context);
    }

    return buffer && buffer.value();
  }

  renderBounds(context) {
    const buffer = context == null ? context = new _path.default() : undefined;
    context.rect(this.xmin, this.ymin, this.xmax - this.xmin, this.ymax - this.ymin);
    return buffer && buffer.value();
  }

  renderCell(i, context) {
    const buffer = context == null ? context = new _path.default() : undefined;

    const points = this._clip(i);

    if (points === null) return;
    context.moveTo(points[0], points[1]);

    for (let i = 2, n = points.length; i < n; i += 2) {
      if (points[i] !== points[i - 2] || points[i + 1] !== points[i - 1]) context.lineTo(points[i], points[i + 1]);
    }

    context.closePath();
    return buffer && buffer.value();
  }

  *cellPolygons() {
    const {
      delaunay: {
        points
      }
    } = this;

    for (let i = 0, n = points.length / 2; i < n; ++i) {
      const cell = this.cellPolygon(i);
      if (cell) yield cell;
    }
  }

  cellPolygon(i) {
    const polygon = new _polygon.default();
    this.renderCell(i, polygon);
    return polygon.value();
  }

  _renderSegment(x0, y0, x1, y1, context) {
    let S;

    const c0 = this._regioncode(x0, y0);

    const c1 = this._regioncode(x1, y1);

    if (c0 === 0 && c1 === 0) {
      context.moveTo(x0, y0);
      context.lineTo(x1, y1);
    } else if (S = this._clipSegment(x0, y0, x1, y1, c0, c1)) {
      context.moveTo(S[0], S[1]);
      context.lineTo(S[2], S[3]);
    }
  }

  contains(i, x, y) {
    if ((x = +x, x !== x) || (y = +y, y !== y)) return false;
    return this.delaunay._step(i, x, y) === i;
  }

  _cell(i) {
    const {
      circumcenters,
      delaunay: {
        inedges,
        halfedges,
        triangles
      }
    } = this;
    const e0 = inedges[i];
    if (e0 === -1) return null; // coincident point

    const points = [];
    let e = e0;

    do {
      const t = Math.floor(e / 3);
      points.push(circumcenters[t * 2], circumcenters[t * 2 + 1]);
      e = e % 3 === 2 ? e - 2 : e + 1;
      if (triangles[e] !== i) break; // bad triangulation

      e = halfedges[e];
    } while (e !== e0 && e !== -1);

    return points;
  }

  _clip(i) {
    // degenerate case (1 valid point: return the box)
    if (i === 0 && this.delaunay.hull.length === 1) {
      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
    }

    const points = this._cell(i);

    if (points === null) return null;
    const {
      vectors: V
    } = this;
    const v = i * 4;
    return V[v] || V[v + 1] ? this._clipInfinite(i, points, V[v], V[v + 1], V[v + 2], V[v + 3]) : this._clipFinite(i, points);
  }

  _clipFinite(i, points) {
    const n = points.length;
    let P = null;
    let x0,
        y0,
        x1 = points[n - 2],
        y1 = points[n - 1];

    let c0,
        c1 = this._regioncode(x1, y1);

    let e0, e1;

    for (let j = 0; j < n; j += 2) {
      x0 = x1, y0 = y1, x1 = points[j], y1 = points[j + 1];
      c0 = c1, c1 = this._regioncode(x1, y1);

      if (c0 === 0 && c1 === 0) {
        e0 = e1, e1 = 0;
        if (P) P.push(x1, y1);else P = [x1, y1];
      } else {
        let S, sx0, sy0, sx1, sy1;

        if (c0 === 0) {
          if ((S = this._clipSegment(x0, y0, x1, y1, c0, c1)) === null) continue;
          [sx0, sy0, sx1, sy1] = S;
        } else {
          if ((S = this._clipSegment(x1, y1, x0, y0, c1, c0)) === null) continue;
          [sx1, sy1, sx0, sy0] = S;
          e0 = e1, e1 = this._edgecode(sx0, sy0);
          if (e0 && e1) this._edge(i, e0, e1, P, P.length);
          if (P) P.push(sx0, sy0);else P = [sx0, sy0];
        }

        e0 = e1, e1 = this._edgecode(sx1, sy1);
        if (e0 && e1) this._edge(i, e0, e1, P, P.length);
        if (P) P.push(sx1, sy1);else P = [sx1, sy1];
      }
    }

    if (P) {
      e0 = e1, e1 = this._edgecode(P[0], P[1]);
      if (e0 && e1) this._edge(i, e0, e1, P, P.length);
    } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
      return [this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax, this.xmin, this.ymin];
    }

    return P;
  }

  _clipSegment(x0, y0, x1, y1, c0, c1) {
    while (true) {
      if (c0 === 0 && c1 === 0) return [x0, y0, x1, y1];
      if (c0 & c1) return null;
      let x,
          y,
          c = c0 || c1;
      if (c & 0b1000) x = x0 + (x1 - x0) * (this.ymax - y0) / (y1 - y0), y = this.ymax;else if (c & 0b0100) x = x0 + (x1 - x0) * (this.ymin - y0) / (y1 - y0), y = this.ymin;else if (c & 0b0010) y = y0 + (y1 - y0) * (this.xmax - x0) / (x1 - x0), x = this.xmax;else y = y0 + (y1 - y0) * (this.xmin - x0) / (x1 - x0), x = this.xmin;
      if (c0) x0 = x, y0 = y, c0 = this._regioncode(x0, y0);else x1 = x, y1 = y, c1 = this._regioncode(x1, y1);
    }
  }

  _clipInfinite(i, points, vx0, vy0, vxn, vyn) {
    let P = Array.from(points),
        p;
    if (p = this._project(P[0], P[1], vx0, vy0)) P.unshift(p[0], p[1]);
    if (p = this._project(P[P.length - 2], P[P.length - 1], vxn, vyn)) P.push(p[0], p[1]);

    if (P = this._clipFinite(i, P)) {
      for (let j = 0, n = P.length, c0, c1 = this._edgecode(P[n - 2], P[n - 1]); j < n; j += 2) {
        c0 = c1, c1 = this._edgecode(P[j], P[j + 1]);
        if (c0 && c1) j = this._edge(i, c0, c1, P, j), n = P.length;
      }
    } else if (this.contains(i, (this.xmin + this.xmax) / 2, (this.ymin + this.ymax) / 2)) {
      P = [this.xmin, this.ymin, this.xmax, this.ymin, this.xmax, this.ymax, this.xmin, this.ymax];
    }

    return P;
  }

  _edge(i, e0, e1, P, j) {
    while (e0 !== e1) {
      let x, y;

      switch (e0) {
        case 0b0101:
          e0 = 0b0100;
          continue;
        // top-left

        case 0b0100:
          e0 = 0b0110, x = this.xmax, y = this.ymin;
          break;
        // top

        case 0b0110:
          e0 = 0b0010;
          continue;
        // top-right

        case 0b0010:
          e0 = 0b1010, x = this.xmax, y = this.ymax;
          break;
        // right

        case 0b1010:
          e0 = 0b1000;
          continue;
        // bottom-right

        case 0b1000:
          e0 = 0b1001, x = this.xmin, y = this.ymax;
          break;
        // bottom

        case 0b1001:
          e0 = 0b0001;
          continue;
        // bottom-left

        case 0b0001:
          e0 = 0b0101, x = this.xmin, y = this.ymin;
          break;
        // left
      }

      if ((P[j] !== x || P[j + 1] !== y) && this.contains(i, x, y)) {
        P.splice(j, 0, x, y), j += 2;
      }
    }

    return j;
  }

  _project(x0, y0, vx, vy) {
    let t = Infinity,
        c,
        x,
        y;

    if (vy < 0) {
      // top
      if (y0 <= this.ymin) return null;
      if ((c = (this.ymin - y0) / vy) < t) y = this.ymin, x = x0 + (t = c) * vx;
    } else if (vy > 0) {
      // bottom
      if (y0 >= this.ymax) return null;
      if ((c = (this.ymax - y0) / vy) < t) y = this.ymax, x = x0 + (t = c) * vx;
    }

    if (vx > 0) {
      // right
      if (x0 >= this.xmax) return null;
      if ((c = (this.xmax - x0) / vx) < t) x = this.xmax, y = y0 + (t = c) * vy;
    } else if (vx < 0) {
      // left
      if (x0 <= this.xmin) return null;
      if ((c = (this.xmin - x0) / vx) < t) x = this.xmin, y = y0 + (t = c) * vy;
    }

    return [x, y];
  }

  _edgecode(x, y) {
    return (x === this.xmin ? 0b0001 : x === this.xmax ? 0b0010 : 0b0000) | (y === this.ymin ? 0b0100 : y === this.ymax ? 0b1000 : 0b0000);
  }

  _regioncode(x, y) {
    return (x < this.xmin ? 0b0001 : x > this.xmax ? 0b0010 : 0b0000) | (y < this.ymin ? 0b0100 : y > this.ymax ? 0b1000 : 0b0000);
  }

}

exports.default = Voronoi;
},{"./path.js":"../node_modules/d3-delaunay/src/path.js","./polygon.js":"../node_modules/d3-delaunay/src/polygon.js"}],"../node_modules/d3-delaunay/src/delaunay.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = void 0;

var _delaunator = _interopRequireDefault(require("delaunator"));

var _path = _interopRequireDefault(require("./path.js"));

var _polygon = _interopRequireDefault(require("./polygon.js"));

var _voronoi = _interopRequireDefault(require("./voronoi.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

const tau = 2 * Math.PI;

function pointX(p) {
  return p[0];
}

function pointY(p) {
  return p[1];
}

function area(hull, points) {
  let n = hull.length,
      x0,
      y0,
      x1 = points[2 * hull[n - 1]],
      y1 = points[2 * hull[n - 1] + 1],
      area = 0;

  for (let i = 0; i < n; i++) {
    x0 = x1, y0 = y1;
    x1 = points[2 * hull[i]];
    y1 = points[2 * hull[i] + 1];
    area += y0 * x1 - x0 * y1;
  }

  return area / 2;
}

function jitter(x, y, r) {
  return [x + Math.sin(x + y) * r, y + Math.cos(x - y) * r];
}

class Delaunay {
  constructor(points) {
    this._delaunator = new _delaunator.default(points);
    this.inedges = new Int32Array(points.length / 2);
    this._hullIndex = new Int32Array(points.length / 2);
    this.points = this._delaunator.coords;

    this._init();
  }

  update() {
    this._delaunator.update();

    this._init();

    return this;
  }

  _init() {
    const d = this._delaunator,
          points = this.points; // check for collinear

    if (d.hull && d.hull.length > 2 && area(d.hull, points) < 1e-10) {
      this.collinear = Int32Array.from({
        length: points.length / 2
      }, (_, i) => i).sort((i, j) => points[2 * i] - points[2 * j] || points[2 * i + 1] - points[2 * j + 1]); // for exact neighbors

      const e = this.collinear[0],
            f = this.collinear[this.collinear.length - 1],
            bounds = [points[2 * e], points[2 * e + 1], points[2 * f], points[2 * f + 1]],
            r = 1e-8 * Math.sqrt((bounds[3] - bounds[1]) ** 2 + (bounds[2] - bounds[0]) ** 2);

      for (let i = 0, n = points.length / 2; i < n; ++i) {
        const p = jitter(points[2 * i], points[2 * i + 1], r);
        points[2 * i] = p[0];
        points[2 * i + 1] = p[1];
      }

      this._delaunator = new _delaunator.default(points);
    } else {
      delete this.collinear;
    }

    const halfedges = this.halfedges = this._delaunator.halfedges;
    const hull = this.hull = this._delaunator.hull;
    const triangles = this.triangles = this._delaunator.triangles;
    const inedges = this.inedges.fill(-1);

    const hullIndex = this._hullIndex.fill(-1); // Compute an index from each point to an (arbitrary) incoming halfedge
    // Used to give the first neighbor of each point; for this reason,
    // on the hull we give priority to exterior halfedges


    for (let e = 0, n = halfedges.length; e < n; ++e) {
      const p = triangles[e % 3 === 2 ? e - 2 : e + 1];
      if (halfedges[e] === -1 || inedges[p] === -1) inedges[p] = e;
    }

    for (let i = 0, n = hull.length; i < n; ++i) {
      hullIndex[hull[i]] = i;
    } // degenerate case: 1 or 2 (distinct) points


    if (hull.length <= 2 && hull.length > 0) {
      this.triangles = new Int32Array(3).fill(-1);
      this.halfedges = new Int32Array(3).fill(-1);
      this.triangles[0] = hull[0];
      this.triangles[1] = hull[1];
      this.triangles[2] = hull[1];
      inedges[hull[0]] = 1;
      if (hull.length === 2) inedges[hull[1]] = 0;
    }
  }

  voronoi(bounds) {
    return new _voronoi.default(this, bounds);
  }

  *neighbors(i) {
    const {
      inedges,
      hull,
      _hullIndex,
      halfedges,
      triangles
    } = this; // degenerate case with several collinear points

    if (this.collinear) {
      const l = this.collinear.indexOf(i);
      if (l > 0) yield this.collinear[l - 1];
      if (l < this.collinear.length - 1) yield this.collinear[l + 1];
      return;
    }

    const e0 = inedges[i];
    if (e0 === -1) return; // coincident point

    let e = e0,
        p0 = -1;

    do {
      yield p0 = triangles[e];
      e = e % 3 === 2 ? e - 2 : e + 1;
      if (triangles[e] !== i) return; // bad triangulation

      e = halfedges[e];

      if (e === -1) {
        const p = hull[(_hullIndex[i] + 1) % hull.length];
        if (p !== p0) yield p;
        return;
      }
    } while (e !== e0);
  }

  find(x, y, i = 0) {
    if ((x = +x, x !== x) || (y = +y, y !== y)) return -1;
    const i0 = i;
    let c;

    while ((c = this._step(i, x, y)) >= 0 && c !== i && c !== i0) i = c;

    return c;
  }

  _step(i, x, y) {
    const {
      inedges,
      hull,
      _hullIndex,
      halfedges,
      triangles,
      points
    } = this;
    if (inedges[i] === -1 || !points.length) return (i + 1) % (points.length >> 1);
    let c = i;
    let dc = (x - points[i * 2]) ** 2 + (y - points[i * 2 + 1]) ** 2;
    const e0 = inedges[i];
    let e = e0;

    do {
      let t = triangles[e];
      const dt = (x - points[t * 2]) ** 2 + (y - points[t * 2 + 1]) ** 2;
      if (dt < dc) dc = dt, c = t;
      e = e % 3 === 2 ? e - 2 : e + 1;
      if (triangles[e] !== i) break; // bad triangulation

      e = halfedges[e];

      if (e === -1) {
        e = hull[(_hullIndex[i] + 1) % hull.length];

        if (e !== t) {
          if ((x - points[e * 2]) ** 2 + (y - points[e * 2 + 1]) ** 2 < dc) return e;
        }

        break;
      }
    } while (e !== e0);

    return c;
  }

  render(context) {
    const buffer = context == null ? context = new _path.default() : undefined;
    const {
      points,
      halfedges,
      triangles
    } = this;

    for (let i = 0, n = halfedges.length; i < n; ++i) {
      const j = halfedges[i];
      if (j < i) continue;
      const ti = triangles[i] * 2;
      const tj = triangles[j] * 2;
      context.moveTo(points[ti], points[ti + 1]);
      context.lineTo(points[tj], points[tj + 1]);
    }

    this.renderHull(context);
    return buffer && buffer.value();
  }

  renderPoints(context, r = 2) {
    const buffer = context == null ? context = new _path.default() : undefined;
    const {
      points
    } = this;

    for (let i = 0, n = points.length; i < n; i += 2) {
      const x = points[i],
            y = points[i + 1];
      context.moveTo(x + r, y);
      context.arc(x, y, r, 0, tau);
    }

    return buffer && buffer.value();
  }

  renderHull(context) {
    const buffer = context == null ? context = new _path.default() : undefined;
    const {
      hull,
      points
    } = this;
    const h = hull[0] * 2,
          n = hull.length;
    context.moveTo(points[h], points[h + 1]);

    for (let i = 1; i < n; ++i) {
      const h = 2 * hull[i];
      context.lineTo(points[h], points[h + 1]);
    }

    context.closePath();
    return buffer && buffer.value();
  }

  hullPolygon() {
    const polygon = new _polygon.default();
    this.renderHull(polygon);
    return polygon.value();
  }

  renderTriangle(i, context) {
    const buffer = context == null ? context = new _path.default() : undefined;
    const {
      points,
      triangles
    } = this;
    const t0 = triangles[i *= 3] * 2;
    const t1 = triangles[i + 1] * 2;
    const t2 = triangles[i + 2] * 2;
    context.moveTo(points[t0], points[t0 + 1]);
    context.lineTo(points[t1], points[t1 + 1]);
    context.lineTo(points[t2], points[t2 + 1]);
    context.closePath();
    return buffer && buffer.value();
  }

  *trianglePolygons() {
    const {
      triangles
    } = this;

    for (let i = 0, n = triangles.length / 3; i < n; ++i) {
      yield this.trianglePolygon(i);
    }
  }

  trianglePolygon(i) {
    const polygon = new _polygon.default();
    this.renderTriangle(i, polygon);
    return polygon.value();
  }

}

exports.default = Delaunay;

Delaunay.from = function (points, fx = pointX, fy = pointY, that) {
  return new Delaunay("length" in points ? flatArray(points, fx, fy, that) : Float64Array.from(flatIterable(points, fx, fy, that)));
};

function flatArray(points, fx, fy, that) {
  const n = points.length;
  const array = new Float64Array(n * 2);

  for (let i = 0; i < n; ++i) {
    const p = points[i];
    array[i * 2] = fx.call(that, p, i, points);
    array[i * 2 + 1] = fy.call(that, p, i, points);
  }

  return array;
}

function* flatIterable(points, fx, fy, that) {
  let i = 0;

  for (const p of points) {
    yield fx.call(that, p, i, points);
    yield fy.call(that, p, i, points);
    ++i;
  }
}
},{"delaunator":"../node_modules/delaunator/index.js","./path.js":"../node_modules/d3-delaunay/src/path.js","./polygon.js":"../node_modules/d3-delaunay/src/polygon.js","./voronoi.js":"../node_modules/d3-delaunay/src/voronoi.js"}],"../node_modules/d3-delaunay/src/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
Object.defineProperty(exports, "Delaunay", {
  enumerable: true,
  get: function () {
    return _delaunay.default;
  }
});
Object.defineProperty(exports, "Voronoi", {
  enumerable: true,
  get: function () {
    return _voronoi.default;
  }
});

var _delaunay = _interopRequireDefault(require("./delaunay.js"));

var _voronoi = _interopRequireDefault(require("./voronoi.js"));

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }
},{"./delaunay.js":"../node_modules/d3-delaunay/src/delaunay.js","./voronoi.js":"../node_modules/d3-delaunay/src/voronoi.js"}],"../node_modules/@ctrl/tinycolor/dist/es/util.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.bound01 = bound01;
exports.clamp01 = clamp01;
exports.isOnePointZero = isOnePointZero;
exports.isPercentage = isPercentage;
exports.boundAlpha = boundAlpha;
exports.convertToPercentage = convertToPercentage;
exports.pad2 = pad2;

function bound01(n, max) {
  if (isOnePointZero(n)) {
    n = '100%';
  }

  var processPercent = isPercentage(n);
  n = max === 360 ? n : Math.min(max, Math.max(0, parseFloat(n)));

  if (processPercent) {
    n = parseInt(String(n * max), 10) / 100;
  }

  if (Math.abs(n - max) < 0.000001) {
    return 1;
  }

  if (max === 360) {
    n = (n < 0 ? n % max + max : n % max) / parseFloat(String(max));
  } else {
    n = n % max / parseFloat(String(max));
  }

  return n;
}

function clamp01(val) {
  return Math.min(1, Math.max(0, val));
}

function isOnePointZero(n) {
  return typeof n === 'string' && n.includes('.') && parseFloat(n) === 1;
}

function isPercentage(n) {
  return typeof n === 'string' && n.includes('%');
}

function boundAlpha(a) {
  a = parseFloat(a);

  if (isNaN(a) || a < 0 || a > 1) {
    a = 1;
  }

  return a;
}

function convertToPercentage(n) {
  if (n <= 1) {
    return Number(n) * 100 + "%";
  }

  return n;
}

function pad2(c) {
  return c.length === 1 ? '0' + c : String(c);
}
},{}],"../node_modules/@ctrl/tinycolor/dist/es/conversion.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.rgbToRgb = rgbToRgb;
exports.rgbToHsl = rgbToHsl;
exports.hslToRgb = hslToRgb;
exports.rgbToHsv = rgbToHsv;
exports.hsvToRgb = hsvToRgb;
exports.rgbToHex = rgbToHex;
exports.rgbaToHex = rgbaToHex;
exports.rgbaToArgbHex = rgbaToArgbHex;
exports.convertDecimalToHex = convertDecimalToHex;
exports.convertHexToDecimal = convertHexToDecimal;
exports.parseIntFromHex = parseIntFromHex;

var _util = require("./util");

function rgbToRgb(r, g, b) {
  return {
    r: (0, _util.bound01)(r, 255) * 255,
    g: (0, _util.bound01)(g, 255) * 255,
    b: (0, _util.bound01)(b, 255) * 255
  };
}

function rgbToHsl(r, g, b) {
  r = (0, _util.bound01)(r, 255);
  g = (0, _util.bound01)(g, 255);
  b = (0, _util.bound01)(b, 255);
  var max = Math.max(r, g, b);
  var min = Math.min(r, g, b);
  var h = 0;
  var s = 0;
  var l = (max + min) / 2;

  if (max === min) {
    s = 0;
    h = 0;
  } else {
    var d = max - min;
    s = l > 0.5 ? d / (2 - max - min) : d / (max + min);

    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;

      case g:
        h = (b - r) / d + 2;
        break;

      case b:
        h = (r - g) / d + 4;
        break;

      default:
        break;
    }

    h /= 6;
  }

  return {
    h: h,
    s: s,
    l: l
  };
}

function hslToRgb(h, s, l) {
  var r;
  var g;
  var b;
  h = (0, _util.bound01)(h, 360);
  s = (0, _util.bound01)(s, 100);
  l = (0, _util.bound01)(l, 100);

  function hue2rgb(p, q, t) {
    if (t < 0) {
      t += 1;
    }

    if (t > 1) {
      t -= 1;
    }

    if (t < 1 / 6) {
      return p + (q - p) * (6 * t);
    }

    if (t < 1 / 2) {
      return q;
    }

    if (t < 2 / 3) {
      return p + (q - p) * (2 / 3 - t) * 6;
    }

    return p;
  }

  if (s === 0) {
    g = l;
    b = l;
    r = l;
  } else {
    var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    var p = 2 * l - q;
    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }

  return {
    r: r * 255,
    g: g * 255,
    b: b * 255
  };
}

function rgbToHsv(r, g, b) {
  r = (0, _util.bound01)(r, 255);
  g = (0, _util.bound01)(g, 255);
  b = (0, _util.bound01)(b, 255);
  var max = Math.max(r, g, b);
  var min = Math.min(r, g, b);
  var h = 0;
  var v = max;
  var d = max - min;
  var s = max === 0 ? 0 : d / max;

  if (max === min) {
    h = 0;
  } else {
    switch (max) {
      case r:
        h = (g - b) / d + (g < b ? 6 : 0);
        break;

      case g:
        h = (b - r) / d + 2;
        break;

      case b:
        h = (r - g) / d + 4;
        break;

      default:
        break;
    }

    h /= 6;
  }

  return {
    h: h,
    s: s,
    v: v
  };
}

function hsvToRgb(h, s, v) {
  h = (0, _util.bound01)(h, 360) * 6;
  s = (0, _util.bound01)(s, 100);
  v = (0, _util.bound01)(v, 100);
  var i = Math.floor(h);
  var f = h - i;
  var p = v * (1 - s);
  var q = v * (1 - f * s);
  var t = v * (1 - (1 - f) * s);
  var mod = i % 6;
  var r = [v, q, p, p, t, v][mod];
  var g = [t, v, v, q, p, p][mod];
  var b = [p, p, t, v, v, q][mod];
  return {
    r: r * 255,
    g: g * 255,
    b: b * 255
  };
}

function rgbToHex(r, g, b, allow3Char) {
  var hex = [(0, _util.pad2)(Math.round(r).toString(16)), (0, _util.pad2)(Math.round(g).toString(16)), (0, _util.pad2)(Math.round(b).toString(16))];

  if (allow3Char && hex[0].charAt(0) === hex[0].charAt(1) && hex[1].charAt(0) === hex[1].charAt(1) && hex[2].charAt(0) === hex[2].charAt(1)) {
    return hex[0].charAt(0) + hex[1].charAt(0) + hex[2].charAt(0);
  }

  return hex.join('');
}

function rgbaToHex(r, g, b, a, allow4Char) {
  var hex = [(0, _util.pad2)(Math.round(r).toString(16)), (0, _util.pad2)(Math.round(g).toString(16)), (0, _util.pad2)(Math.round(b).toString(16)), (0, _util.pad2)(convertDecimalToHex(a))];

  if (allow4Char && hex[0].charAt(0) === hex[0].charAt(1) && hex[1].charAt(0) === hex[1].charAt(1) && hex[2].charAt(0) === hex[2].charAt(1) && hex[3].charAt(0) === hex[3].charAt(1)) {
    return hex[0].charAt(0) + hex[1].charAt(0) + hex[2].charAt(0) + hex[3].charAt(0);
  }

  return hex.join('');
}

function rgbaToArgbHex(r, g, b, a) {
  var hex = [(0, _util.pad2)(convertDecimalToHex(a)), (0, _util.pad2)(Math.round(r).toString(16)), (0, _util.pad2)(Math.round(g).toString(16)), (0, _util.pad2)(Math.round(b).toString(16))];
  return hex.join('');
}

function convertDecimalToHex(d) {
  return Math.round(parseFloat(d) * 255).toString(16);
}

function convertHexToDecimal(h) {
  return parseIntFromHex(h) / 255;
}

function parseIntFromHex(val) {
  return parseInt(val, 16);
}
},{"./util":"../node_modules/@ctrl/tinycolor/dist/es/util.js"}],"../node_modules/@ctrl/tinycolor/dist/es/css-color-names.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.names = void 0;
var names = {
  aliceblue: '#f0f8ff',
  antiquewhite: '#faebd7',
  aqua: '#00ffff',
  aquamarine: '#7fffd4',
  azure: '#f0ffff',
  beige: '#f5f5dc',
  bisque: '#ffe4c4',
  black: '#000000',
  blanchedalmond: '#ffebcd',
  blue: '#0000ff',
  blueviolet: '#8a2be2',
  brown: '#a52a2a',
  burlywood: '#deb887',
  cadetblue: '#5f9ea0',
  chartreuse: '#7fff00',
  chocolate: '#d2691e',
  coral: '#ff7f50',
  cornflowerblue: '#6495ed',
  cornsilk: '#fff8dc',
  crimson: '#dc143c',
  cyan: '#00ffff',
  darkblue: '#00008b',
  darkcyan: '#008b8b',
  darkgoldenrod: '#b8860b',
  darkgray: '#a9a9a9',
  darkgreen: '#006400',
  darkgrey: '#a9a9a9',
  darkkhaki: '#bdb76b',
  darkmagenta: '#8b008b',
  darkolivegreen: '#556b2f',
  darkorange: '#ff8c00',
  darkorchid: '#9932cc',
  darkred: '#8b0000',
  darksalmon: '#e9967a',
  darkseagreen: '#8fbc8f',
  darkslateblue: '#483d8b',
  darkslategray: '#2f4f4f',
  darkslategrey: '#2f4f4f',
  darkturquoise: '#00ced1',
  darkviolet: '#9400d3',
  deeppink: '#ff1493',
  deepskyblue: '#00bfff',
  dimgray: '#696969',
  dimgrey: '#696969',
  dodgerblue: '#1e90ff',
  firebrick: '#b22222',
  floralwhite: '#fffaf0',
  forestgreen: '#228b22',
  fuchsia: '#ff00ff',
  gainsboro: '#dcdcdc',
  ghostwhite: '#f8f8ff',
  gold: '#ffd700',
  goldenrod: '#daa520',
  gray: '#808080',
  green: '#008000',
  greenyellow: '#adff2f',
  grey: '#808080',
  honeydew: '#f0fff0',
  hotpink: '#ff69b4',
  indianred: '#cd5c5c',
  indigo: '#4b0082',
  ivory: '#fffff0',
  khaki: '#f0e68c',
  lavender: '#e6e6fa',
  lavenderblush: '#fff0f5',
  lawngreen: '#7cfc00',
  lemonchiffon: '#fffacd',
  lightblue: '#add8e6',
  lightcoral: '#f08080',
  lightcyan: '#e0ffff',
  lightgoldenrodyellow: '#fafad2',
  lightgray: '#d3d3d3',
  lightgreen: '#90ee90',
  lightgrey: '#d3d3d3',
  lightpink: '#ffb6c1',
  lightsalmon: '#ffa07a',
  lightseagreen: '#20b2aa',
  lightskyblue: '#87cefa',
  lightslategray: '#778899',
  lightslategrey: '#778899',
  lightsteelblue: '#b0c4de',
  lightyellow: '#ffffe0',
  lime: '#00ff00',
  limegreen: '#32cd32',
  linen: '#faf0e6',
  magenta: '#ff00ff',
  maroon: '#800000',
  mediumaquamarine: '#66cdaa',
  mediumblue: '#0000cd',
  mediumorchid: '#ba55d3',
  mediumpurple: '#9370db',
  mediumseagreen: '#3cb371',
  mediumslateblue: '#7b68ee',
  mediumspringgreen: '#00fa9a',
  mediumturquoise: '#48d1cc',
  mediumvioletred: '#c71585',
  midnightblue: '#191970',
  mintcream: '#f5fffa',
  mistyrose: '#ffe4e1',
  moccasin: '#ffe4b5',
  navajowhite: '#ffdead',
  navy: '#000080',
  oldlace: '#fdf5e6',
  olive: '#808000',
  olivedrab: '#6b8e23',
  orange: '#ffa500',
  orangered: '#ff4500',
  orchid: '#da70d6',
  palegoldenrod: '#eee8aa',
  palegreen: '#98fb98',
  paleturquoise: '#afeeee',
  palevioletred: '#db7093',
  papayawhip: '#ffefd5',
  peachpuff: '#ffdab9',
  peru: '#cd853f',
  pink: '#ffc0cb',
  plum: '#dda0dd',
  powderblue: '#b0e0e6',
  purple: '#800080',
  rebeccapurple: '#663399',
  red: '#ff0000',
  rosybrown: '#bc8f8f',
  royalblue: '#4169e1',
  saddlebrown: '#8b4513',
  salmon: '#fa8072',
  sandybrown: '#f4a460',
  seagreen: '#2e8b57',
  seashell: '#fff5ee',
  sienna: '#a0522d',
  silver: '#c0c0c0',
  skyblue: '#87ceeb',
  slateblue: '#6a5acd',
  slategray: '#708090',
  slategrey: '#708090',
  snow: '#fffafa',
  springgreen: '#00ff7f',
  steelblue: '#4682b4',
  tan: '#d2b48c',
  teal: '#008080',
  thistle: '#d8bfd8',
  tomato: '#ff6347',
  turquoise: '#40e0d0',
  violet: '#ee82ee',
  wheat: '#f5deb3',
  white: '#ffffff',
  whitesmoke: '#f5f5f5',
  yellow: '#ffff00',
  yellowgreen: '#9acd32'
};
exports.names = names;
},{}],"../node_modules/@ctrl/tinycolor/dist/es/format-input.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.inputToRGB = inputToRGB;
exports.stringInputToObject = stringInputToObject;
exports.isValidCSSUnit = isValidCSSUnit;

var _conversion = require("./conversion");

var _cssColorNames = require("./css-color-names");

var _util = require("./util");

function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

function inputToRGB(color) {
  var rgb = {
    r: 0,
    g: 0,
    b: 0
  };
  var a = 1;
  var s = null;
  var v = null;
  var l = null;
  var ok = false;
  var format = false;

  if (typeof color === 'string') {
    color = stringInputToObject(color);
  }

  if (_typeof(color) === 'object') {
    if (isValidCSSUnit(color.r) && isValidCSSUnit(color.g) && isValidCSSUnit(color.b)) {
      rgb = (0, _conversion.rgbToRgb)(color.r, color.g, color.b);
      ok = true;
      format = String(color.r).substr(-1) === '%' ? 'prgb' : 'rgb';
    } else if (isValidCSSUnit(color.h) && isValidCSSUnit(color.s) && isValidCSSUnit(color.v)) {
      s = (0, _util.convertToPercentage)(color.s);
      v = (0, _util.convertToPercentage)(color.v);
      rgb = (0, _conversion.hsvToRgb)(color.h, s, v);
      ok = true;
      format = 'hsv';
    } else if (isValidCSSUnit(color.h) && isValidCSSUnit(color.s) && isValidCSSUnit(color.l)) {
      s = (0, _util.convertToPercentage)(color.s);
      l = (0, _util.convertToPercentage)(color.l);
      rgb = (0, _conversion.hslToRgb)(color.h, s, l);
      ok = true;
      format = 'hsl';
    }

    if (Object.prototype.hasOwnProperty.call(color, 'a')) {
      a = color.a;
    }
  }

  a = (0, _util.boundAlpha)(a);
  return {
    ok: ok,
    format: color.format || format,
    r: Math.min(255, Math.max(rgb.r, 0)),
    g: Math.min(255, Math.max(rgb.g, 0)),
    b: Math.min(255, Math.max(rgb.b, 0)),
    a: a
  };
}

var CSS_INTEGER = '[-\\+]?\\d+%?';
var CSS_NUMBER = '[-\\+]?\\d*\\.\\d+%?';
var CSS_UNIT = "(?:" + CSS_NUMBER + ")|(?:" + CSS_INTEGER + ")";
var PERMISSIVE_MATCH3 = "[\\s|\\(]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")\\s*\\)?";
var PERMISSIVE_MATCH4 = "[\\s|\\(]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")[,|\\s]+(" + CSS_UNIT + ")\\s*\\)?";
var matchers = {
  CSS_UNIT: new RegExp(CSS_UNIT),
  rgb: new RegExp('rgb' + PERMISSIVE_MATCH3),
  rgba: new RegExp('rgba' + PERMISSIVE_MATCH4),
  hsl: new RegExp('hsl' + PERMISSIVE_MATCH3),
  hsla: new RegExp('hsla' + PERMISSIVE_MATCH4),
  hsv: new RegExp('hsv' + PERMISSIVE_MATCH3),
  hsva: new RegExp('hsva' + PERMISSIVE_MATCH4),
  hex3: /^#?([0-9a-fA-F]{1})([0-9a-fA-F]{1})([0-9a-fA-F]{1})$/,
  hex6: /^#?([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})$/,
  hex4: /^#?([0-9a-fA-F]{1})([0-9a-fA-F]{1})([0-9a-fA-F]{1})([0-9a-fA-F]{1})$/,
  hex8: /^#?([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})([0-9a-fA-F]{2})$/
};

function stringInputToObject(color) {
  color = color.trim().toLowerCase();

  if (color.length === 0) {
    return false;
  }

  var named = false;

  if (_cssColorNames.names[color]) {
    color = _cssColorNames.names[color];
    named = true;
  } else if (color === 'transparent') {
    return {
      r: 0,
      g: 0,
      b: 0,
      a: 0,
      format: 'name'
    };
  }

  var match = matchers.rgb.exec(color);

  if (match) {
    return {
      r: match[1],
      g: match[2],
      b: match[3]
    };
  }

  match = matchers.rgba.exec(color);

  if (match) {
    return {
      r: match[1],
      g: match[2],
      b: match[3],
      a: match[4]
    };
  }

  match = matchers.hsl.exec(color);

  if (match) {
    return {
      h: match[1],
      s: match[2],
      l: match[3]
    };
  }

  match = matchers.hsla.exec(color);

  if (match) {
    return {
      h: match[1],
      s: match[2],
      l: match[3],
      a: match[4]
    };
  }

  match = matchers.hsv.exec(color);

  if (match) {
    return {
      h: match[1],
      s: match[2],
      v: match[3]
    };
  }

  match = matchers.hsva.exec(color);

  if (match) {
    return {
      h: match[1],
      s: match[2],
      v: match[3],
      a: match[4]
    };
  }

  match = matchers.hex8.exec(color);

  if (match) {
    return {
      r: (0, _conversion.parseIntFromHex)(match[1]),
      g: (0, _conversion.parseIntFromHex)(match[2]),
      b: (0, _conversion.parseIntFromHex)(match[3]),
      a: (0, _conversion.convertHexToDecimal)(match[4]),
      format: named ? 'name' : 'hex8'
    };
  }

  match = matchers.hex6.exec(color);

  if (match) {
    return {
      r: (0, _conversion.parseIntFromHex)(match[1]),
      g: (0, _conversion.parseIntFromHex)(match[2]),
      b: (0, _conversion.parseIntFromHex)(match[3]),
      format: named ? 'name' : 'hex'
    };
  }

  match = matchers.hex4.exec(color);

  if (match) {
    return {
      r: (0, _conversion.parseIntFromHex)(match[1] + match[1]),
      g: (0, _conversion.parseIntFromHex)(match[2] + match[2]),
      b: (0, _conversion.parseIntFromHex)(match[3] + match[3]),
      a: (0, _conversion.convertHexToDecimal)(match[4] + match[4]),
      format: named ? 'name' : 'hex8'
    };
  }

  match = matchers.hex3.exec(color);

  if (match) {
    return {
      r: (0, _conversion.parseIntFromHex)(match[1] + match[1]),
      g: (0, _conversion.parseIntFromHex)(match[2] + match[2]),
      b: (0, _conversion.parseIntFromHex)(match[3] + match[3]),
      format: named ? 'name' : 'hex'
    };
  }

  return false;
}

function isValidCSSUnit(color) {
  return Boolean(matchers.CSS_UNIT.exec(String(color)));
}
},{"./conversion":"../node_modules/@ctrl/tinycolor/dist/es/conversion.js","./css-color-names":"../node_modules/@ctrl/tinycolor/dist/es/css-color-names.js","./util":"../node_modules/@ctrl/tinycolor/dist/es/util.js"}],"../node_modules/@ctrl/tinycolor/dist/es/index.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.tinycolor = tinycolor;
exports.TinyColor = void 0;

var _conversion = require("./conversion");

var _cssColorNames = require("./css-color-names");

var _formatInput = require("./format-input");

var _util = require("./util");

var TinyColor = function () {
  function TinyColor(color, opts) {
    if (color === void 0) {
      color = '';
    }

    if (opts === void 0) {
      opts = {};
    }

    if (color instanceof TinyColor) {
      return color;
    }

    this.originalInput = color;
    var rgb = (0, _formatInput.inputToRGB)(color);
    this.originalInput = color;
    this.r = rgb.r;
    this.g = rgb.g;
    this.b = rgb.b;
    this.a = rgb.a;
    this.roundA = Math.round(100 * this.a) / 100;
    this.format = opts.format || rgb.format;
    this.gradientType = opts.gradientType;

    if (this.r < 1) {
      this.r = Math.round(this.r);
    }

    if (this.g < 1) {
      this.g = Math.round(this.g);
    }

    if (this.b < 1) {
      this.b = Math.round(this.b);
    }

    this.isValid = rgb.ok;
  }

  TinyColor.prototype.isDark = function () {
    return this.getBrightness() < 128;
  };

  TinyColor.prototype.isLight = function () {
    return !this.isDark();
  };

  TinyColor.prototype.getBrightness = function () {
    var rgb = this.toRgb();
    return (rgb.r * 299 + rgb.g * 587 + rgb.b * 114) / 1000;
  };

  TinyColor.prototype.getLuminance = function () {
    var rgb = this.toRgb();
    var R;
    var G;
    var B;
    var RsRGB = rgb.r / 255;
    var GsRGB = rgb.g / 255;
    var BsRGB = rgb.b / 255;

    if (RsRGB <= 0.03928) {
      R = RsRGB / 12.92;
    } else {
      R = Math.pow((RsRGB + 0.055) / 1.055, 2.4);
    }

    if (GsRGB <= 0.03928) {
      G = GsRGB / 12.92;
    } else {
      G = Math.pow((GsRGB + 0.055) / 1.055, 2.4);
    }

    if (BsRGB <= 0.03928) {
      B = BsRGB / 12.92;
    } else {
      B = Math.pow((BsRGB + 0.055) / 1.055, 2.4);
    }

    return 0.2126 * R + 0.7152 * G + 0.0722 * B;
  };

  TinyColor.prototype.getAlpha = function () {
    return this.a;
  };

  TinyColor.prototype.setAlpha = function (alpha) {
    this.a = (0, _util.boundAlpha)(alpha);
    this.roundA = Math.round(100 * this.a) / 100;
    return this;
  };

  TinyColor.prototype.toHsv = function () {
    var hsv = (0, _conversion.rgbToHsv)(this.r, this.g, this.b);
    return {
      h: hsv.h * 360,
      s: hsv.s,
      v: hsv.v,
      a: this.a
    };
  };

  TinyColor.prototype.toHsvString = function () {
    var hsv = (0, _conversion.rgbToHsv)(this.r, this.g, this.b);
    var h = Math.round(hsv.h * 360);
    var s = Math.round(hsv.s * 100);
    var v = Math.round(hsv.v * 100);
    return this.a === 1 ? "hsv(" + h + ", " + s + "%, " + v + "%)" : "hsva(" + h + ", " + s + "%, " + v + "%, " + this.roundA + ")";
  };

  TinyColor.prototype.toHsl = function () {
    var hsl = (0, _conversion.rgbToHsl)(this.r, this.g, this.b);
    return {
      h: hsl.h * 360,
      s: hsl.s,
      l: hsl.l,
      a: this.a
    };
  };

  TinyColor.prototype.toHslString = function () {
    var hsl = (0, _conversion.rgbToHsl)(this.r, this.g, this.b);
    var h = Math.round(hsl.h * 360);
    var s = Math.round(hsl.s * 100);
    var l = Math.round(hsl.l * 100);
    return this.a === 1 ? "hsl(" + h + ", " + s + "%, " + l + "%)" : "hsla(" + h + ", " + s + "%, " + l + "%, " + this.roundA + ")";
  };

  TinyColor.prototype.toHex = function (allow3Char) {
    if (allow3Char === void 0) {
      allow3Char = false;
    }

    return (0, _conversion.rgbToHex)(this.r, this.g, this.b, allow3Char);
  };

  TinyColor.prototype.toHexString = function (allow3Char) {
    if (allow3Char === void 0) {
      allow3Char = false;
    }

    return '#' + this.toHex(allow3Char);
  };

  TinyColor.prototype.toHex8 = function (allow4Char) {
    if (allow4Char === void 0) {
      allow4Char = false;
    }

    return (0, _conversion.rgbaToHex)(this.r, this.g, this.b, this.a, allow4Char);
  };

  TinyColor.prototype.toHex8String = function (allow4Char) {
    if (allow4Char === void 0) {
      allow4Char = false;
    }

    return '#' + this.toHex8(allow4Char);
  };

  TinyColor.prototype.toRgb = function () {
    return {
      r: Math.round(this.r),
      g: Math.round(this.g),
      b: Math.round(this.b),
      a: this.a
    };
  };

  TinyColor.prototype.toRgbString = function () {
    var r = Math.round(this.r);
    var g = Math.round(this.g);
    var b = Math.round(this.b);
    return this.a === 1 ? "rgb(" + r + ", " + g + ", " + b + ")" : "rgba(" + r + ", " + g + ", " + b + ", " + this.roundA + ")";
  };

  TinyColor.prototype.toPercentageRgb = function () {
    var fmt = function fmt(x) {
      return Math.round((0, _util.bound01)(x, 255) * 100) + "%";
    };

    return {
      r: fmt(this.r),
      g: fmt(this.g),
      b: fmt(this.b),
      a: this.a
    };
  };

  TinyColor.prototype.toPercentageRgbString = function () {
    var rnd = function rnd(x) {
      return Math.round((0, _util.bound01)(x, 255) * 100);
    };

    return this.a === 1 ? "rgb(" + rnd(this.r) + "%, " + rnd(this.g) + "%, " + rnd(this.b) + "%)" : "rgba(" + rnd(this.r) + "%, " + rnd(this.g) + "%, " + rnd(this.b) + "%, " + this.roundA + ")";
  };

  TinyColor.prototype.toName = function () {
    if (this.a === 0) {
      return 'transparent';
    }

    if (this.a < 1) {
      return false;
    }

    var hex = '#' + (0, _conversion.rgbToHex)(this.r, this.g, this.b, false);

    for (var _i = 0, _a = Object.keys(_cssColorNames.names); _i < _a.length; _i++) {
      var key = _a[_i];

      if (_cssColorNames.names[key] === hex) {
        return key;
      }
    }

    return false;
  };

  TinyColor.prototype.toString = function (format) {
    var formatSet = Boolean(format);
    format = format || this.format;
    var formattedString = false;
    var hasAlpha = this.a < 1 && this.a >= 0;
    var needsAlphaFormat = !formatSet && hasAlpha && (format.startsWith('hex') || format === 'name');

    if (needsAlphaFormat) {
      if (format === 'name' && this.a === 0) {
        return this.toName();
      }

      return this.toRgbString();
    }

    if (format === 'rgb') {
      formattedString = this.toRgbString();
    }

    if (format === 'prgb') {
      formattedString = this.toPercentageRgbString();
    }

    if (format === 'hex' || format === 'hex6') {
      formattedString = this.toHexString();
    }

    if (format === 'hex3') {
      formattedString = this.toHexString(true);
    }

    if (format === 'hex4') {
      formattedString = this.toHex8String(true);
    }

    if (format === 'hex8') {
      formattedString = this.toHex8String();
    }

    if (format === 'name') {
      formattedString = this.toName();
    }

    if (format === 'hsl') {
      formattedString = this.toHslString();
    }

    if (format === 'hsv') {
      formattedString = this.toHsvString();
    }

    return formattedString || this.toHexString();
  };

  TinyColor.prototype.clone = function () {
    return new TinyColor(this.toString());
  };

  TinyColor.prototype.lighten = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    var hsl = this.toHsl();
    hsl.l += amount / 100;
    hsl.l = (0, _util.clamp01)(hsl.l);
    return new TinyColor(hsl);
  };

  TinyColor.prototype.brighten = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    var rgb = this.toRgb();
    rgb.r = Math.max(0, Math.min(255, rgb.r - Math.round(255 * -(amount / 100))));
    rgb.g = Math.max(0, Math.min(255, rgb.g - Math.round(255 * -(amount / 100))));
    rgb.b = Math.max(0, Math.min(255, rgb.b - Math.round(255 * -(amount / 100))));
    return new TinyColor(rgb);
  };

  TinyColor.prototype.darken = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    var hsl = this.toHsl();
    hsl.l -= amount / 100;
    hsl.l = (0, _util.clamp01)(hsl.l);
    return new TinyColor(hsl);
  };

  TinyColor.prototype.tint = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    return this.mix('white', amount);
  };

  TinyColor.prototype.shade = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    return this.mix('black', amount);
  };

  TinyColor.prototype.desaturate = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    var hsl = this.toHsl();
    hsl.s -= amount / 100;
    hsl.s = (0, _util.clamp01)(hsl.s);
    return new TinyColor(hsl);
  };

  TinyColor.prototype.saturate = function (amount) {
    if (amount === void 0) {
      amount = 10;
    }

    var hsl = this.toHsl();
    hsl.s += amount / 100;
    hsl.s = (0, _util.clamp01)(hsl.s);
    return new TinyColor(hsl);
  };

  TinyColor.prototype.greyscale = function () {
    return this.desaturate(100);
  };

  TinyColor.prototype.spin = function (amount) {
    var hsl = this.toHsl();
    var hue = (hsl.h + amount) % 360;
    hsl.h = hue < 0 ? 360 + hue : hue;
    return new TinyColor(hsl);
  };

  TinyColor.prototype.mix = function (color, amount) {
    if (amount === void 0) {
      amount = 50;
    }

    var rgb1 = this.toRgb();
    var rgb2 = new TinyColor(color).toRgb();
    var p = amount / 100;
    var rgba = {
      r: (rgb2.r - rgb1.r) * p + rgb1.r,
      g: (rgb2.g - rgb1.g) * p + rgb1.g,
      b: (rgb2.b - rgb1.b) * p + rgb1.b,
      a: (rgb2.a - rgb1.a) * p + rgb1.a
    };
    return new TinyColor(rgba);
  };

  TinyColor.prototype.analogous = function (results, slices) {
    if (results === void 0) {
      results = 6;
    }

    if (slices === void 0) {
      slices = 30;
    }

    var hsl = this.toHsl();
    var part = 360 / slices;
    var ret = [this];

    for (hsl.h = (hsl.h - (part * results >> 1) + 720) % 360; --results;) {
      hsl.h = (hsl.h + part) % 360;
      ret.push(new TinyColor(hsl));
    }

    return ret;
  };

  TinyColor.prototype.complement = function () {
    var hsl = this.toHsl();
    hsl.h = (hsl.h + 180) % 360;
    return new TinyColor(hsl);
  };

  TinyColor.prototype.monochromatic = function (results) {
    if (results === void 0) {
      results = 6;
    }

    var hsv = this.toHsv();
    var h = hsv.h;
    var s = hsv.s;
    var v = hsv.v;
    var res = [];
    var modification = 1 / results;

    while (results--) {
      res.push(new TinyColor({
        h: h,
        s: s,
        v: v
      }));
      v = (v + modification) % 1;
    }

    return res;
  };

  TinyColor.prototype.splitcomplement = function () {
    var hsl = this.toHsl();
    var h = hsl.h;
    return [this, new TinyColor({
      h: (h + 72) % 360,
      s: hsl.s,
      l: hsl.l
    }), new TinyColor({
      h: (h + 216) % 360,
      s: hsl.s,
      l: hsl.l
    })];
  };

  TinyColor.prototype.triad = function () {
    return this.polyad(3);
  };

  TinyColor.prototype.tetrad = function () {
    return this.polyad(4);
  };

  TinyColor.prototype.polyad = function (n) {
    var hsl = this.toHsl();
    var h = hsl.h;
    var result = [this];
    var increment = 360 / n;

    for (var i = 1; i < n; i++) {
      result.push(new TinyColor({
        h: (h + i * increment) % 360,
        s: hsl.s,
        l: hsl.l
      }));
    }

    return result;
  };

  TinyColor.prototype.equals = function (color) {
    return this.toRgbString() === new TinyColor(color).toRgbString();
  };

  return TinyColor;
}();

exports.TinyColor = TinyColor;

function tinycolor(color, opts) {
  if (color === void 0) {
    color = '';
  }

  if (opts === void 0) {
    opts = {};
  }

  return new TinyColor(color, opts);
}
},{"./conversion":"../node_modules/@ctrl/tinycolor/dist/es/conversion.js","./css-color-names":"../node_modules/@ctrl/tinycolor/dist/es/css-color-names.js","./format-input":"../node_modules/@ctrl/tinycolor/dist/es/format-input.js","./util":"../node_modules/@ctrl/tinycolor/dist/es/util.js"}],"../node_modules/@ctrl/tinycolor/dist/es/readability.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.readability = readability;
exports.isReadable = isReadable;
exports.mostReadable = mostReadable;

var _index = require("./index");

function readability(color1, color2) {
  var c1 = new _index.TinyColor(color1);
  var c2 = new _index.TinyColor(color2);
  return (Math.max(c1.getLuminance(), c2.getLuminance()) + 0.05) / (Math.min(c1.getLuminance(), c2.getLuminance()) + 0.05);
}

function isReadable(color1, color2, wcag2) {
  if (wcag2 === void 0) {
    wcag2 = {
      level: 'AA',
      size: 'small'
    };
  }

  var readabilityLevel = readability(color1, color2);

  switch ((wcag2.level || 'AA') + (wcag2.size || 'small')) {
    case 'AAsmall':
    case 'AAAlarge':
      return readabilityLevel >= 4.5;

    case 'AAlarge':
      return readabilityLevel >= 3;

    case 'AAAsmall':
      return readabilityLevel >= 7;

    default:
      return false;
  }
}

function mostReadable(baseColor, colorList, args) {
  if (args === void 0) {
    args = {
      includeFallbackColors: false,
      level: 'AA',
      size: 'small'
    };
  }

  var bestColor = null;
  var bestScore = 0;
  var includeFallbackColors = args.includeFallbackColors,
      level = args.level,
      size = args.size;

  for (var _i = 0, colorList_1 = colorList; _i < colorList_1.length; _i++) {
    var color = colorList_1[_i];
    var score = readability(baseColor, color);

    if (score > bestScore) {
      bestScore = score;
      bestColor = new _index.TinyColor(color);
    }
  }

  if (isReadable(baseColor, bestColor, {
    level: level,
    size: size
  }) || !includeFallbackColors) {
    return bestColor;
  }

  args.includeFallbackColors = false;
  return mostReadable(baseColor, ['#fff', '#000'], args);
}
},{"./index":"../node_modules/@ctrl/tinycolor/dist/es/index.js"}],"../node_modules/@ctrl/tinycolor/dist/es/to-ms-filter.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.toMsFilter = toMsFilter;

var _conversion = require("./conversion");

var _index = require("./index");

function toMsFilter(firstColor, secondColor) {
  var color = new _index.TinyColor(firstColor);
  var hex8String = '#' + (0, _conversion.rgbaToArgbHex)(color.r, color.g, color.b, color.a);
  var secondHex8String = hex8String;
  var gradientType = color.gradientType ? 'GradientType = 1, ' : '';

  if (secondColor) {
    var s = new _index.TinyColor(secondColor);
    secondHex8String = '#' + (0, _conversion.rgbaToArgbHex)(s.r, s.g, s.b, s.a);
  }

  return "progid:DXImageTransform.Microsoft.gradient(" + gradientType + "startColorstr=" + hex8String + ",endColorstr=" + secondHex8String + ")";
}
},{"./conversion":"../node_modules/@ctrl/tinycolor/dist/es/conversion.js","./index":"../node_modules/@ctrl/tinycolor/dist/es/index.js"}],"../node_modules/@ctrl/tinycolor/dist/es/from-ratio.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.fromRatio = fromRatio;
exports.legacyRandom = legacyRandom;

var _index = require("./index");

var _util = require("./util");

function fromRatio(ratio, opts) {
  var newColor = {
    r: (0, _util.convertToPercentage)(ratio.r),
    g: (0, _util.convertToPercentage)(ratio.g),
    b: (0, _util.convertToPercentage)(ratio.b)
  };

  if (ratio.a !== undefined) {
    newColor.a = Number(ratio.a);
  }

  return new _index.TinyColor(newColor, opts);
}

function legacyRandom() {
  return new _index.TinyColor({
    r: Math.random(),
    g: Math.random(),
    b: Math.random()
  });
}
},{"./index":"../node_modules/@ctrl/tinycolor/dist/es/index.js","./util":"../node_modules/@ctrl/tinycolor/dist/es/util.js"}],"../node_modules/@ctrl/tinycolor/dist/es/random.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.random = random;
exports.bounds = void 0;

var _index = require("./index");

function random(options) {
  if (options === void 0) {
    options = {};
  }

  if (options.count !== undefined && options.count !== null) {
    var totalColors = options.count;
    var colors = [];
    options.count = undefined;

    while (totalColors > colors.length) {
      options.count = null;

      if (options.seed) {
        options.seed += 1;
      }

      colors.push(random(options));
    }

    options.count = totalColors;
    return colors;
  }

  var h = pickHue(options.hue, options.seed);
  var s = pickSaturation(h, options);
  var v = pickBrightness(h, s, options);
  var res = {
    h: h,
    s: s,
    v: v
  };

  if (options.alpha !== undefined) {
    res.a = options.alpha;
  }

  return new _index.TinyColor(res);
}

function pickHue(hue, seed) {
  var hueRange = getHueRange(hue);
  var res = randomWithin(hueRange, seed);

  if (res < 0) {
    res = 360 + res;
  }

  return res;
}

function pickSaturation(hue, options) {
  if (options.hue === 'monochrome') {
    return 0;
  }

  if (options.luminosity === 'random') {
    return randomWithin([0, 100], options.seed);
  }

  var saturationRange = getColorInfo(hue).saturationRange;
  var sMin = saturationRange[0];
  var sMax = saturationRange[1];

  switch (options.luminosity) {
    case 'bright':
      sMin = 55;
      break;

    case 'dark':
      sMin = sMax - 10;
      break;

    case 'light':
      sMax = 55;
      break;

    default:
      break;
  }

  return randomWithin([sMin, sMax], options.seed);
}

function pickBrightness(H, S, options) {
  var bMin = getMinimumBrightness(H, S);
  var bMax = 100;

  switch (options.luminosity) {
    case 'dark':
      bMax = bMin + 20;
      break;

    case 'light':
      bMin = (bMax + bMin) / 2;
      break;

    case 'random':
      bMin = 0;
      bMax = 100;
      break;

    default:
      break;
  }

  return randomWithin([bMin, bMax], options.seed);
}

function getMinimumBrightness(H, S) {
  var lowerBounds = getColorInfo(H).lowerBounds;

  for (var i = 0; i < lowerBounds.length - 1; i++) {
    var s1 = lowerBounds[i][0];
    var v1 = lowerBounds[i][1];
    var s2 = lowerBounds[i + 1][0];
    var v2 = lowerBounds[i + 1][1];

    if (S >= s1 && S <= s2) {
      var m = (v2 - v1) / (s2 - s1);
      var b = v1 - m * s1;
      return m * S + b;
    }
  }

  return 0;
}

function getHueRange(colorInput) {
  var num = parseInt(colorInput, 10);

  if (!Number.isNaN(num) && num < 360 && num > 0) {
    return [num, num];
  }

  if (typeof colorInput === 'string') {
    var namedColor = bounds.find(function (n) {
      return n.name === colorInput;
    });

    if (namedColor) {
      var color = defineColor(namedColor);

      if (color.hueRange) {
        return color.hueRange;
      }
    }

    var parsed = new _index.TinyColor(colorInput);

    if (parsed.isValid) {
      var hue = parsed.toHsv().h;
      return [hue, hue];
    }
  }

  return [0, 360];
}

function getColorInfo(hue) {
  if (hue >= 334 && hue <= 360) {
    hue -= 360;
  }

  for (var _i = 0, bounds_1 = bounds; _i < bounds_1.length; _i++) {
    var bound = bounds_1[_i];
    var color = defineColor(bound);

    if (color.hueRange && hue >= color.hueRange[0] && hue <= color.hueRange[1]) {
      return color;
    }
  }

  throw Error('Color not found');
}

function randomWithin(range, seed) {
  if (seed === undefined) {
    return Math.floor(range[0] + Math.random() * (range[1] + 1 - range[0]));
  }

  var max = range[1] || 1;
  var min = range[0] || 0;
  seed = (seed * 9301 + 49297) % 233280;
  var rnd = seed / 233280.0;
  return Math.floor(min + rnd * (max - min));
}

function defineColor(bound) {
  var sMin = bound.lowerBounds[0][0];
  var sMax = bound.lowerBounds[bound.lowerBounds.length - 1][0];
  var bMin = bound.lowerBounds[bound.lowerBounds.length - 1][1];
  var bMax = bound.lowerBounds[0][1];
  return {
    name: bound.name,
    hueRange: bound.hueRange,
    lowerBounds: bound.lowerBounds,
    saturationRange: [sMin, sMax],
    brightnessRange: [bMin, bMax]
  };
}

var bounds = [{
  name: 'monochrome',
  hueRange: null,
  lowerBounds: [[0, 0], [100, 0]]
}, {
  name: 'red',
  hueRange: [-26, 18],
  lowerBounds: [[20, 100], [30, 92], [40, 89], [50, 85], [60, 78], [70, 70], [80, 60], [90, 55], [100, 50]]
}, {
  name: 'orange',
  hueRange: [19, 46],
  lowerBounds: [[20, 100], [30, 93], [40, 88], [50, 86], [60, 85], [70, 70], [100, 70]]
}, {
  name: 'yellow',
  hueRange: [47, 62],
  lowerBounds: [[25, 100], [40, 94], [50, 89], [60, 86], [70, 84], [80, 82], [90, 80], [100, 75]]
}, {
  name: 'green',
  hueRange: [63, 178],
  lowerBounds: [[30, 100], [40, 90], [50, 85], [60, 81], [70, 74], [80, 64], [90, 50], [100, 40]]
}, {
  name: 'blue',
  hueRange: [179, 257],
  lowerBounds: [[20, 100], [30, 86], [40, 80], [50, 74], [60, 60], [70, 52], [80, 44], [90, 39], [100, 35]]
}, {
  name: 'purple',
  hueRange: [258, 282],
  lowerBounds: [[20, 100], [30, 87], [40, 79], [50, 70], [60, 65], [70, 59], [80, 52], [90, 45], [100, 42]]
}, {
  name: 'pink',
  hueRange: [283, 334],
  lowerBounds: [[20, 100], [30, 90], [40, 86], [60, 84], [80, 80], [90, 75], [100, 73]]
}];
exports.bounds = bounds;
},{"./index":"../node_modules/@ctrl/tinycolor/dist/es/index.js"}],"../node_modules/@ctrl/tinycolor/dist/es/public_api.js":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
var _exportNames = {};
exports.default = void 0;

var _index = require("./index");

Object.keys(_index).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _index[key];
    }
  });
});

var _cssColorNames = require("./css-color-names");

Object.keys(_cssColorNames).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _cssColorNames[key];
    }
  });
});

var _readability = require("./readability");

Object.keys(_readability).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _readability[key];
    }
  });
});

var _toMsFilter = require("./to-ms-filter");

Object.keys(_toMsFilter).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _toMsFilter[key];
    }
  });
});

var _fromRatio = require("./from-ratio");

Object.keys(_fromRatio).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _fromRatio[key];
    }
  });
});

var _formatInput = require("./format-input");

Object.keys(_formatInput).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _formatInput[key];
    }
  });
});

var _random = require("./random");

Object.keys(_random).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function () {
      return _random[key];
    }
  });
});
var _default = _index.tinycolor;
exports.default = _default;
},{"./index":"../node_modules/@ctrl/tinycolor/dist/es/index.js","./css-color-names":"../node_modules/@ctrl/tinycolor/dist/es/css-color-names.js","./readability":"../node_modules/@ctrl/tinycolor/dist/es/readability.js","./to-ms-filter":"../node_modules/@ctrl/tinycolor/dist/es/to-ms-filter.js","./from-ratio":"../node_modules/@ctrl/tinycolor/dist/es/from-ratio.js","./format-input":"../node_modules/@ctrl/tinycolor/dist/es/format-input.js","./random":"../node_modules/@ctrl/tinycolor/dist/es/random.js"}],"main.ts":[function(require,module,exports) {
"use strict";

exports.__esModule = true;

var vec2_1 = require("./math/vec2");

var d3_delaunay_1 = require("d3-delaunay");

var tinycolor_1 = require("@ctrl/tinycolor");

var NUM_POINTS = 5000;
var WIDTH = 1280;
var HEIGHT = 720;
var NUM_PLATES = 10;
var PLATE_COLOR_MAP = ["black", "blue", "red", "purple", "green", "teal", "navy", "yellow", "magenta", "white", "pink"];
var NUM_PLATE_STEPS = 10;
var PLATE_ANGLE_COS = 0;
var MAX_HEIGHT_DIFF = 5;
var LINELESS = true;
var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");

function choose(arr) {
  return Math.floor(Math.random() * arr.length);
}

function generatePoints() {
  var res = new Array(NUM_POINTS);

  for (var i = 0; i < NUM_POINTS; i++) {
    res[i] = {
      loc: {
        x: Math.random() * WIDTH,
        y: Math.random() * HEIGHT
      },
      plateId: -1,
      height: 0
    };
  }

  return res;
}

function dist(p1, p2) {
  return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
}

function floodFill(points) {
  var unmapped_points = points.slice();
  var res = [];

  for (var i = 0; i < NUM_PLATES; i++) {
    var p_i = choose(unmapped_points);
    var p = unmapped_points.splice(p_i, 1)[0];
    p.plateId = i;
    res.push(p);
  }

  while (unmapped_points.length > 0) {
    var p_i = choose(unmapped_points);
    var p = unmapped_points.splice(p_i, 1)[0];
    var closest_i = -1;
    var closest_dist = 0;
    var rl = res.length;

    for (var i = 0; i < rl; i++) {
      var curr = res[i];
      var d = dist(curr.loc, p.loc);

      if (closest_i == -1 || d < closest_dist) {
        closest_dist = d;
        closest_i = i;
      }
    }

    p.plateId = res[closest_i].plateId;
    res.push(p);
  }

  return res;
}

function genNeighbors(voronoi) {
  var res = [];
  var cells = [];

  for (var i = 0; i < NUM_POINTS; i++) {
    cells.push(voronoi.cellPolygon(i));
  }

  var _loop_1 = function _loop_1(i) {
    var cell_poly = voronoi.cellPolygon(i);
    var neighbor_cells = cells.map(function (val, idx) {
      return {
        val: val,
        idx: idx
      };
    }).filter(function (_a) {
      var val = _a.val;
      var matched_verts = val.filter(function (pnt) {
        return cell_poly.filter(function (pnt2) {
          return pnt[0] == pnt2[0] && pnt[1] == pnt[1];
        }).length > 0;
      }).length;
      return matched_verts == 2;
    }).map(function (val) {
      return val.idx;
    });
    res.push(neighbor_cells);
  };

  for (var i = 0; i < NUM_POINTS; i++) {
    _loop_1(i);
  }

  return res;
}

function getRandom(start, end) {
  return Math.floor(Math.random() * (end - start + 1) + start);
}

function shuffleInPlace(array) {
  var _a; // if it's 1 or 0 items, just return


  if (array.length <= 1) return array; // For each index in array

  for (var i = 0; i < array.length; i++) {
    // choose a random not-yet-placed item to place there
    // must be an item AFTER the current item, because the stuff
    // before has all already been placed
    var randomChoiceIndex = getRandom(i, array.length - 1); // place our random choice in the spot by swapping

    _a = [array[randomChoiceIndex], array[i]], array[i] = _a[0], array[randomChoiceIndex] = _a[1];
  }

  return array;
}

function range(start, end) {
  return Array.from({
    length: end - start
  }, function (v, k) {
    return k + start;
  });
}

function heightGen(voronoi, points) {
  var plateDirections = [];

  for (var plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
    plateDirections.push(new vec2_1.Vec2(Math.random() * 2 - 1, Math.random() * 2 - 1));
  }

  var neighbors = genNeighbors(voronoi);
  var r = range(0, NUM_POINTS);

  for (var i = 0; i < NUM_PLATE_STEPS; i++) {
    shuffleInPlace(r);

    for (var _i = 0, r_1 = r; _i < r_1.length; _i++) {
      var cell = r_1[_i];
      var cell_point = points[cell];
      var neighbor_cells = neighbors[cell];
      shuffleInPlace(neighbor_cells);

      for (var _a = 0, neighbor_cells_1 = neighbor_cells; _a < neighbor_cells_1.length; _a++) {
        var neighbor_cell_idx = neighbor_cells_1[_a];
        var neighbor_cell_point = points[neighbor_cell_idx];
        var dir = new vec2_1.Vec2(neighbor_cell_point.loc.x - cell_point.loc.x, neighbor_cell_point.loc.y - cell_point.loc.y);
        var cell_plate_dir = plateDirections[cell_point.plateId];
        var dot = dir.x * cell_plate_dir.x + dir.y * cell_plate_dir.y;

        if (dot > PLATE_ANGLE_COS) {
          neighbor_cell_point.height++;
          cell_point.height++;
        } else if (dot < -PLATE_ANGLE_COS) {
          neighbor_cell_point.height--;
          cell_point.height--;
        }

        if (Math.abs(cell_point.height - neighbor_cell_point.height) > MAX_HEIGHT_DIFF) {
          var avg = (cell_point.height + neighbor_cell_point.height) / 2;
          cell_point.height = (cell_point.height + avg) / 2;
          neighbor_cell_point.height = (neighbor_cell_point.height + avg) / 2;
        }
      }
    }
  }

  return points;
}

function main() {
  var points = generatePoints();
  floodFill(points);
  var voronoi_points = points.map(function (val) {
    return [val.loc.x, val.loc.y];
  });
  var delaunay = d3_delaunay_1.Delaunay.from(voronoi_points);
  var voronoi = delaunay.voronoi([0, 0, WIDTH, HEIGHT]);

  for (var i = 0; i < NUM_POINTS; i++) {
    var cp = voronoi.cellPolygon(i);
    points[i].loc.x = cp.reduce(function (agg, val) {
      return agg + val[0];
    }, 0) / cp.length;
    points[i].loc.y = cp.reduce(function (agg, val) {
      return agg + val[1];
    }, 0) / cp.length;
  }

  voronoi_points = points.map(function (val) {
    return [val.loc.x, val.loc.y];
  });
  delaunay = d3_delaunay_1.Delaunay.from(voronoi_points);
  voronoi = delaunay.voronoi([0, 0, WIDTH, HEIGHT]);
  heightGen(voronoi, points);
  var min_height = points.reduce(function (agg, curr) {
    return curr.height < agg ? curr.height : agg;
  }, points[0].height);
  var max_height = points.reduce(function (agg, curr) {
    return curr.height > agg ? curr.height : agg;
  }, points[0].height);
  var avg_height = points.reduce(function (agg, curr) {
    return agg + curr.height;
  }, 0) / points.length;
  ctx.clearRect(0, 0, 800, 600);
  ctx.strokeStyle = "black";

  var _loop_2 = function _loop_2(plate_idx) {
    console.log("Rendering " + ctx.fillStyle);
    var pts = points.map(function (val, idx) {
      return {
        p: val,
        idx: idx
      };
    }).filter(function (val) {
      return val.p.plateId === plate_idx;
    });
    var pl = pts.length;

    for (var i = 0; i < pl; i++) {
      ctx.beginPath();
      var p = pts[i];
      voronoi.renderCell(p.idx, ctx);
      var height_percent = 100 - 100 * (p.p.height - min_height) / (max_height - min_height);
      var color = void 0;

      if (p.p.height > avg_height) {
        color = new tinycolor_1.TinyColor('brown');
      } else {
        color = new tinycolor_1.TinyColor('blue');
      }

      color = color.shade(height_percent / 1.2);
      ctx.fillStyle = color.toRgbString(); // ctx.fillStyle = `rgb(${height_percent}%, ${height_percent}%, ${height_percent}%)`

      ctx.fill();

      if (LINELESS) {
        ctx.strokeStyle = ctx.fillStyle;
      }

      ctx.stroke();
    }
  };

  for (var plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
    _loop_2(plate_idx);
  }
}

main();
},{"./math/vec2":"math/vec2.ts","d3-delaunay":"../node_modules/d3-delaunay/src/index.js","@ctrl/tinycolor":"../node_modules/@ctrl/tinycolor/dist/es/public_api.js"}],"../node_modules/parcel-bundler/src/builtins/hmr-runtime.js":[function(require,module,exports) {
var global = arguments[3];
var OVERLAY_ID = '__parcel__error__overlay__';
var OldModule = module.bundle.Module;

function Module(moduleName) {
  OldModule.call(this, moduleName);
  this.hot = {
    data: module.bundle.hotData,
    _acceptCallbacks: [],
    _disposeCallbacks: [],
    accept: function (fn) {
      this._acceptCallbacks.push(fn || function () {});
    },
    dispose: function (fn) {
      this._disposeCallbacks.push(fn);
    }
  };
  module.bundle.hotData = null;
}

module.bundle.Module = Module;
var checkedAssets, assetsToAccept;
var parent = module.bundle.parent;

if ((!parent || !parent.isParcelRequire) && typeof WebSocket !== 'undefined') {
  var hostname = "" || location.hostname;
  var protocol = location.protocol === 'https:' ? 'wss' : 'ws';
  var ws = new WebSocket(protocol + '://' + hostname + ':' + "50947" + '/');

  ws.onmessage = function (event) {
    checkedAssets = {};
    assetsToAccept = [];
    var data = JSON.parse(event.data);

    if (data.type === 'update') {
      var handled = false;
      data.assets.forEach(function (asset) {
        if (!asset.isNew) {
          var didAccept = hmrAcceptCheck(global.parcelRequire, asset.id);

          if (didAccept) {
            handled = true;
          }
        }
      }); // Enable HMR for CSS by default.

      handled = handled || data.assets.every(function (asset) {
        return asset.type === 'css' && asset.generated.js;
      });

      if (handled) {
        console.clear();
        data.assets.forEach(function (asset) {
          hmrApply(global.parcelRequire, asset);
        });
        assetsToAccept.forEach(function (v) {
          hmrAcceptRun(v[0], v[1]);
        });
      } else {
        window.location.reload();
      }
    }

    if (data.type === 'reload') {
      ws.close();

      ws.onclose = function () {
        location.reload();
      };
    }

    if (data.type === 'error-resolved') {
      console.log('[parcel] ✨ Error resolved');
      removeErrorOverlay();
    }

    if (data.type === 'error') {
      console.error('[parcel] 🚨  ' + data.error.message + '\n' + data.error.stack);
      removeErrorOverlay();
      var overlay = createErrorOverlay(data);
      document.body.appendChild(overlay);
    }
  };
}

function removeErrorOverlay() {
  var overlay = document.getElementById(OVERLAY_ID);

  if (overlay) {
    overlay.remove();
  }
}

function createErrorOverlay(data) {
  var overlay = document.createElement('div');
  overlay.id = OVERLAY_ID; // html encode message and stack trace

  var message = document.createElement('div');
  var stackTrace = document.createElement('pre');
  message.innerText = data.error.message;
  stackTrace.innerText = data.error.stack;
  overlay.innerHTML = '<div style="background: black; font-size: 16px; color: white; position: fixed; height: 100%; width: 100%; top: 0px; left: 0px; padding: 30px; opacity: 0.85; font-family: Menlo, Consolas, monospace; z-index: 9999;">' + '<span style="background: red; padding: 2px 4px; border-radius: 2px;">ERROR</span>' + '<span style="top: 2px; margin-left: 5px; position: relative;">🚨</span>' + '<div style="font-size: 18px; font-weight: bold; margin-top: 20px;">' + message.innerHTML + '</div>' + '<pre>' + stackTrace.innerHTML + '</pre>' + '</div>';
  return overlay;
}

function getParents(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return [];
  }

  var parents = [];
  var k, d, dep;

  for (k in modules) {
    for (d in modules[k][1]) {
      dep = modules[k][1][d];

      if (dep === id || Array.isArray(dep) && dep[dep.length - 1] === id) {
        parents.push(k);
      }
    }
  }

  if (bundle.parent) {
    parents = parents.concat(getParents(bundle.parent, id));
  }

  return parents;
}

function hmrApply(bundle, asset) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (modules[asset.id] || !bundle.parent) {
    var fn = new Function('require', 'module', 'exports', asset.generated.js);
    asset.isNew = !modules[asset.id];
    modules[asset.id] = [fn, asset.deps];
  } else if (bundle.parent) {
    hmrApply(bundle.parent, asset);
  }
}

function hmrAcceptCheck(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (!modules[id] && bundle.parent) {
    return hmrAcceptCheck(bundle.parent, id);
  }

  if (checkedAssets[id]) {
    return;
  }

  checkedAssets[id] = true;
  var cached = bundle.cache[id];
  assetsToAccept.push([bundle, id]);

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    return true;
  }

  return getParents(global.parcelRequire, id).some(function (id) {
    return hmrAcceptCheck(global.parcelRequire, id);
  });
}

function hmrAcceptRun(bundle, id) {
  var cached = bundle.cache[id];
  bundle.hotData = {};

  if (cached) {
    cached.hot.data = bundle.hotData;
  }

  if (cached && cached.hot && cached.hot._disposeCallbacks.length) {
    cached.hot._disposeCallbacks.forEach(function (cb) {
      cb(bundle.hotData);
    });
  }

  delete bundle.cache[id];
  bundle(id);
  cached = bundle.cache[id];

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    cached.hot._acceptCallbacks.forEach(function (cb) {
      cb();
    });

    return true;
  }
}
},{}]},{},["../node_modules/parcel-bundler/src/builtins/hmr-runtime.js","main.ts"], null)
//# sourceMappingURL=/main.c39d6dcf.js.map