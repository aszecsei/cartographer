export function lerp(a0: number, a1: number, t: number): number {
    return (1.0 - t)*a0 + t*a1
}