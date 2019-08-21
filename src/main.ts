import { Vec2 } from "./math/vec2";
import { Delaunay, Voronoi } from "d3-delaunay";
import { lerp } from './math/mathf';

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
const NUM_PLATE_STEPS = 10;

interface IPoint {
  loc: Vec2;
  plateId: number;
  height: number;
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
      plateId: -1,
      height: 0,
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

function genNeighbors(voronoi: Voronoi<number[]>): number[][] {
  let res = [];

  let cells: number[][][] = [];
  for (let i = 0; i < NUM_POINTS; i++) {
    cells.push(voronoi.cellPolygon(i));
  }

  for (let i = 0; i < NUM_POINTS; i++) {
    let cell_poly = voronoi.cellPolygon(i);
    let neighbor_cells = cells.map((val, idx) => ({val, idx})).filter(({val}) => {
      let matched_verts = val.filter(pnt => cell_poly.filter(pnt2 => pnt[0] == pnt2[0] && pnt[1] == pnt[1]).length > 0).length;
      return matched_verts == 2
    }).map(val => val.idx);
    res.push(neighbor_cells);
  }
  
  return res;
}

function heightGen(voronoi: Voronoi<number[]>, points: IPoint[]): IPoint[] {
  let plateDirections: Vec2[] = [];
  for (let plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
    plateDirections.push(new Vec2(
      Math.random() * 2 - 1,
       Math.random() * 2 - 1
    ));
  }

  let neighbors = genNeighbors(voronoi);

  for (let i = 0; i < NUM_PLATE_STEPS; i++) {
    for (let cell = 0; cell < NUM_POINTS; cell++) {
      let cell_point = points[cell];
      let neighbor_cells = neighbors[cell];
      for (var neighbor_cell_idx of neighbor_cells) {
        let neighbor_cell_point = points[neighbor_cell_idx];
        let dir = new Vec2(
          neighbor_cell_point.loc.x - cell_point.loc.x,
          neighbor_cell_point.loc.y - cell_point.loc.y,
        );
        let cell_plate_dir = plateDirections[cell_point.plateId];
        let dot = dir.x * cell_plate_dir.x + dir.y * cell_plate_dir.y;
        if (dot > 0) {
          if (neighbor_cell_point.height - cell_point.height <= 3) {
            neighbor_cell_point.height++;
          } else {
            cell_point.height++;
          }
        } else {
          if (cell_point.height - neighbor_cell_point.height <= 3) {
            neighbor_cell_point.height--;
          } else {
            cell_point.height--;
          }
        }
      }
    }
  }

  return points;
}

function main() {
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

  heightGen(voronoi, points);

  let min_height = points.reduce((agg, curr) => (curr.height < agg ? curr.height : agg), points[0].height);
  let max_height = points.reduce((agg, curr) => (curr.height > agg ? curr.height : agg), points[0].height);

  ctx.clearRect(0, 0, 800, 600);

  ctx.strokeStyle = "black";
  for (let plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
    console.log("Rendering " + ctx.fillStyle);
    let pts = points
      .map((val, idx) => ({ p: val, idx }))
      .filter(val => val.p.plateId === plate_idx);
    let pl = pts.length;
    for (let i = 0; i < pl; i++) {
      ctx.beginPath();
      
      let p = pts[i];
      voronoi.renderCell(p.idx, ctx);

      let height_percent = (p.p.height - min_height) / (max_height - min_height) * 100;
      ctx.fillStyle = `rgb(${height_percent}%, ${height_percent}%, ${height_percent}%)`
      ctx.fill();
      ctx.stroke();
    }

    // ctx.fillStyle = PLATE_COLOR_MAP[plate_idx];
    // ctx.fill();
    // ctx.stroke();
    
  }
}

main();