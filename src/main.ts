import { Vec2 } from "./math/vec2";
import { Delaunay, Voronoi } from "d3-delaunay";

const NUM_POINTS = 5000;
const WIDTH = 1280;
const HEIGHT = 720;

const NUM_PLATES = 10;
const PLATE_COLOR_MAP = [
  "black",
  "blue",
  "red",
  "purple",
  "green",
  "teal",
  "navy",
  "yellow",
  "magenta",
  "white",
  "pink"
];

interface IPoint {
  loc: Vec2;
  plateId: number;
}

let canvas = document.getElementById("canvas") as HTMLCanvasElement;
let ctx = canvas.getContext("2d");

function choose<T>(arr: Array<T>): number {
  return Math.floor(Math.random() * arr.length);
}

function generatePoints(): IPoint[] {
  let res = new Array<IPoint>(NUM_POINTS);
  for (let i = 0; i < NUM_POINTS; i++) {
    res[i] = {
      loc: {
        x: Math.random() * WIDTH,
        y: Math.random() * HEIGHT
      },
      plateId: -1
    };
  }
  return res;
}

function dist(p1: Vec2, p2: Vec2): number {
  return Math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2);
}

function floodFill(points: IPoint[]): IPoint[] {
  let unmapped_points = [...points];
  let res = [];
  for (let i = 0; i < NUM_PLATES; i++) {
    let p_i = choose(unmapped_points);
    let p = unmapped_points.splice(p_i, 1)[0];
    p.plateId = i;
    res.push(p);
  }
  while (unmapped_points.length > 0) {
    let p_i = choose(unmapped_points);
    let p = unmapped_points.splice(p_i, 1)[0];
    let closest_i = -1;
    let closest_dist = 0;
    let rl = res.length;
    for (let i = 0; i < rl; i++) {
      let curr = res[i];
      let d = dist(curr.loc, p.loc);
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

let points = generatePoints();
floodFill(points);
let voronoi_points = points.map(val => [val.loc.x, val.loc.y]);
let delaunay = Delaunay.from(voronoi_points);
let voronoi = delaunay.voronoi([0, 0, WIDTH, HEIGHT]);

for (let i = 0; i < NUM_POINTS; i++) {
  let cp = voronoi.cellPolygon(i);
  points[i].loc.x = cp.reduce((agg, val) => agg + val[0], 0) / cp.length;
  points[i].loc.y = cp.reduce((agg, val) => agg + val[1], 0) / cp.length;
}

voronoi_points = points.map(val => [val.loc.x, val.loc.y]);
delaunay = Delaunay.from(voronoi_points);
voronoi = delaunay.voronoi([0, 0, WIDTH, HEIGHT]);

ctx.clearRect(0, 0, 800, 600);

ctx.strokeStyle = "black";
for (let plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
  ctx.beginPath();
  console.log("Rendering " + ctx.fillStyle);
  let pts = points
    .map((val, idx) => ({ p: val, idx }))
    .filter(val => val.p.plateId === plate_idx);
  let pl = pts.length;
  for (let i = 0; i < pl; i++) {
    let p = pts[i];
    voronoi.renderCell(p.idx, ctx);
  }
  ctx.fillStyle = PLATE_COLOR_MAP[plate_idx];
  ctx.fill();
  ctx.stroke();
}
