import { Vec2 } from "./math/vec2";
import { Delaunay, Voronoi } from "d3-delaunay";
import { lerp } from './math/mathf';
import { TinyColor } from '@ctrl/tinycolor';

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
const NUM_PLATE_STEPS = 20;
const PLATE_ANGLE_COS = 0;
const MAX_HEIGHT_DIFF = 5;

const WATER_ITERATIONS = 10000;
const RAIN_TRANSFER_AMT = 1;
const EROSION_AMT = 0.2;
const CYCLES_PER_EVAPORATION = 1000;
const EVAPORATION_AMT = 0.2;

const LINELESS = true;

interface IPoint {
  loc: Vec2;
  plateId: number;
  height: number;
  water: number;
}

function getTotalHeight(point: IPoint): number {
  return point.height + point.water;
}

function smoothStep(a: number, b: number, t: number): number {
  let new_t = t * t * (3 - 2 * t);
  return lerp(a, b, new_t);
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
      water: 0,
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
    let cell_poly = cells[i];
    let neighbor_cells = cells.map((val, idx) => ({val, idx})).filter(({val, idx}) => {
      let matched_verts = val.filter(pnt => cell_poly.filter(pnt2 => pnt[0] == pnt2[0] && pnt[1] == pnt[1]).length > 0).length;
      return (idx != i) && (matched_verts >= 2);
    }).map(val => val.idx);
    if (neighbor_cells.length == 0) {
      console.error("Cell has no neighbors!");
    }
    res.push(neighbor_cells);
  }
  
  return res;
}

function getRandom(start, end) {
  return Math.floor(Math.random() * (end - start + 1) + start)
}

function shuffleInPlace<T>(array: T[]): T[] {
  // if it's 1 or 0 items, just return
  if (array.length <= 1) return array;

  // For each index in array
  for (let i = 0; i < array.length; i++) {

    // choose a random not-yet-placed item to place there
    // must be an item AFTER the current item, because the stuff
    // before has all already been placed
    const randomChoiceIndex = getRandom(i, array.length - 1);

    // place our random choice in the spot by swapping
    [array[i], array[randomChoiceIndex]] = [array[randomChoiceIndex], array[i]];
  }

  return array;
}


function range(start, end): number[] { return Array.from({length: (end - start)}, (v, k) => k + start); }

function heightGen(voronoi: Voronoi<number[]>, points: IPoint[]): IPoint[] {
  let plateDirections: Vec2[] = [];
  for (let plate_idx = 0; plate_idx < NUM_PLATES; plate_idx++) {
    plateDirections.push(new Vec2(
      Math.random() * 2 - 1,
       Math.random() * 2 - 1
    ));
  }

  let neighbors = genNeighbors(voronoi);

  let r = range(0, NUM_POINTS);

  for (let i = 0; i < NUM_PLATE_STEPS; i++) {
    shuffleInPlace(r);
    for (let cell of r) {
      let cell_point = points[cell];
      let neighbor_cells = [...neighbors[cell]];
      shuffleInPlace(neighbor_cells);
      for (var neighbor_cell_idx of neighbor_cells) {
        let neighbor_cell_point = points[neighbor_cell_idx];
        let dir = new Vec2(
          neighbor_cell_point.loc.x - cell_point.loc.x,
          neighbor_cell_point.loc.y - cell_point.loc.y,
        );
        let cell_plate_dir = plateDirections[cell_point.plateId];
        let dot = dir.x * cell_plate_dir.x + dir.y * cell_plate_dir.y;
        if (dot > PLATE_ANGLE_COS) {
          neighbor_cell_point.height++;
          cell_point.height++;
        } else if (dot < -PLATE_ANGLE_COS) {
          neighbor_cell_point.height--;
          cell_point.height--;
        }

        if (Math.abs(cell_point.height - neighbor_cell_point.height) > MAX_HEIGHT_DIFF) {
          let avg = (cell_point.height + neighbor_cell_point.height) / 2;
          cell_point.height = (cell_point.height + avg) / 2;
          neighbor_cell_point.height = (neighbor_cell_point.height + avg) / 2;
        }
      }
    }
  }

  return points;
}

function waterGen(voronoi: Voronoi<number[]>, points: IPoint[]): IPoint[] {
  let neighbors = genNeighbors(voronoi);

  let avg_height = points.reduce((agg, curr) => (curr.height + agg), points[0].height) / NUM_POINTS;
  for (let r of range(0, NUM_POINTS)) {
    let cell = points[r];
    if (cell.height < avg_height) {
      cell.water = avg_height - cell.height;
    }
  }

  let r = range(0, NUM_POINTS);

  for (let i = 0; i < WATER_ITERATIONS; i++) {
    shuffleInPlace(r);

    // Evaporation (only occurs every CYCLES_PER_EVAPORATION iteration)
    if (i % CYCLES_PER_EVAPORATION == 0) {
      let total_evaporation = 0;

      for (let cell of r) {
        let cell_point = points[cell];
        if (cell_point.water > 0) {
          let amt = (cell_point.water > EVAPORATION_AMT ? EVAPORATION_AMT : cell_point.water);
          cell_point.water -= amt;
          total_evaporation += amt;
        }
      }

      let rainfall = total_evaporation / NUM_POINTS;
      console.log("RAINFALL: " + rainfall);

      for (let cell of r) {
        let cell_point = points[cell];
        cell_point.water += rainfall;
      }
    }

    for (let cell of r) {
      let cell_point = points[cell];
      let neighbor_cells = [...neighbors[cell]];
      shuffleInPlace(neighbor_cells);
            
      if (cell_point.water > 0) {
        // Water flows from a cell to its lowest neighbor if the cell
        // (a) has water, and
        // (b) the cell's height plus water is greater than the neighbor cell's height plus water
        let min_neighbor = neighbor_cells.reduce((curr_min, neighbor_idx) => {
          if (curr_min == undefined) { return points[neighbor_idx]; }
          if (getTotalHeight(points[neighbor_idx]) < getTotalHeight(curr_min)) { return points[neighbor_idx] }
          else { return curr_min }
        }, undefined as IPoint | undefined);
        if (min_neighbor === undefined) {
          console.error("NEIGHBOR UNDEFINED???");
          console.log(neighbor_cells);
          continue;
        }

        let transfer_diff = getTotalHeight(cell_point) - getTotalHeight(min_neighbor);
        if (transfer_diff > 0) {
          if (transfer_diff <= RAIN_TRANSFER_AMT) {
            // Average the two heights
            let avg_height = ((cell_point.height + cell_point.water) + (min_neighbor.height + min_neighbor.water)) / 2;
            let transfer_amt = avg_height - (min_neighbor.height + min_neighbor.water);

            if (transfer_amt < cell_point.water) {
              transfer_amt = cell_point.water;
            }

            cell_point.water -= transfer_amt;
            min_neighbor.water += transfer_amt;

            let percent = transfer_amt / RAIN_TRANSFER_AMT;
            // Erosion
            cell_point.height -= EROSION_AMT * transfer_diff * percent;
            min_neighbor.height += EROSION_AMT * transfer_diff * percent;
          } else {
            cell_point.water -= RAIN_TRANSFER_AMT;
            min_neighbor.water += RAIN_TRANSFER_AMT;

            // Erosion
            cell_point.height -= EROSION_AMT * transfer_diff;
            min_neighbor.height += EROSION_AMT * transfer_diff;
          }
        }
      }
    }
  }

  return points;
}

function render(voronoi: Voronoi<number[]>, points: IPoint[]) {
  let min_height = points.reduce((agg, curr) => (getTotalHeight(curr) < agg ? getTotalHeight(curr) : agg), getTotalHeight(points[0]));
  let max_height = points.reduce((agg, curr) => (getTotalHeight(curr) > agg ? getTotalHeight(curr) : agg), getTotalHeight(points[0]));

  let max_water = points.reduce((agg, curr) => (curr.water > agg ? curr.water : agg), points[0].water);

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

      let height_percent = 100 - (100 * (getTotalHeight(p.p) - min_height) / (max_height - min_height));
      let color: TinyColor;
      if (p.p.water < EVAPORATION_AMT) {
        color = new TinyColor('brown');
        color = color.shade(height_percent / 1.1);
      } else {
        let water_percent = lerp(0, 100, (p.p.water / max_water));
        let shallow = new TinyColor("#0059ff");
        let deep = new TinyColor("#2109ab");
        color = shallow.mix(deep, water_percent);
        // color = new TinyColor("blue");
      }
      
      ctx.fillStyle = color.toRgbString();
      // ctx.fillStyle = `rgb(${height_percent}%, ${height_percent}%, ${height_percent}%)`
      ctx.fill();

      if (LINELESS) { ctx.strokeStyle = ctx.fillStyle; }
      ctx.stroke();
    }

    // ctx.fillStyle = PLATE_COLOR_MAP[plate_idx];
    // ctx.fill();
    // ctx.stroke();
    
  }
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
  waterGen(voronoi, points);

  render(voronoi, points);
}

main();