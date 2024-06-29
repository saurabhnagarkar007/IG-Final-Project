//ensures requestAnimationFrame is available in all browsers. 
//If not available, it falls back to setTimeout to simulate a frame rate of approximately 60 frames per second, crucial for smooth animation.
window.requestAnimFrame = window.requestAnimationFrame ||
function (callback) {
  window.setTimeout(callback, 1e3 / 120);
};

let gravity = 100; //Downward force applied to each point.

// Dimensions of the cloth, defining how many points are created vertically and horizontally.
let clothY = 24; 
let clothX = 70; 

let spacing = 18; // Distance between points.
let tearDist = 60;// Maximum length a link can stretch before it breaks.

// Control damping and how points react upon hitting the canvas edges.
let friction = 0.99;
let bounce = 0.5;

// Canvas setup.
let canvas = document.getElementById('canvas');
let ctx = canvas.getContext('2d');
canvas.width = Math.min(1920, window.innerWidth);
canvas.height = 750;
ctx.strokeStyle = 'blue';

// Mouse setup.
let mouse = {
  cut: 8,
  influence: 36,
  down: false,
  button: 1,
  x: 0,
  y: 0,
  px: 0,
  py: 0 };

//Point class models a single point in the cloth. Each point can interact with physical forces, connect with other points, and be pinned in space
class Point {
  constructor(x, y) {
    this.x = x;
    this.y = y;
    this.px = x;
    this.py = y;
    this.vx = 0;
    this.vy = 0;
    this.pinX = null;
    this.pinY = null;

    this.constraints = [];
  }

  update(delta) { //Calculates new position using Verlet integration, considering gravity and constraints.
    if (this.pinX && this.pinY) return this;

    // Handling mouse influence or cutting
    if (mouse.down) {
      let dx = this.x - mouse.x;   // Distance between mouse x and point x
      let dy = this.y - mouse.y;   // Distance between mouse y and point y
      let dist = Math.sqrt(dx * dx + dy * dy);  // Calculating Euclidean distance

      // Mouse dragging influence
      if (mouse.button === 1 && dist < mouse.influence) {
        this.px = this.x - (mouse.x - mouse.px);
        this.py = this.y - (mouse.y - mouse.py);

      } else if (dist < mouse.cut) {     // Mouse cutting action
        this.constraints = [];    // Removes all constraints, effectively cutting the cloth
      }
    }

    // Applying gravity force only in the y direction
    this.addForce(0, gravity);

    // Verlet Integration for position update
    let nx = this.x + (this.x - this.px) * friction + this.vx * delta;
    let ny = this.y + (this.y - this.py) * friction + this.vy * delta;

    // Updating velocities to zero after applying them (damping)
    this.px = this.x;
    this.py = this.y;

    this.x = nx;
    this.y = ny;

    this.vy = this.vx = 0;

    // Boundary conditions checking and bouncing logic
    if (this.x >= canvas.width) {
      this.px = canvas.width + (canvas.width - this.px) * bounce;
      this.x = canvas.width;
    } else if (this.x <= 0) {
      this.px *= -1 * bounce;
      this.x = 0;
    }

    if (this.y >= canvas.height) {
      this.py = canvas.height + (canvas.height - this.py) * bounce;
      this.y = canvas.height;
    } else if (this.y <= 0) {
      this.py *= -1 * bounce;
      this.y = 0;
    }

    return this;
  }

  draw() {   //Draws lines to other points as defined in its constraints.
    let i = this.constraints.length;
    while (i--) this.constraints[i].draw();
  }

  resolve() {     //Iterates over each constraint and resolves them to maintain the structure of the cloth.
    if (this.pinX && this.pinY) {
      this.x = this.pinX;
      this.y = this.pinY;
      return;
    }

    this.constraints.forEach(constraint => constraint.resolve());
  }

  attach(point) {   //Adds a new constraint between this point and another.
    this.constraints.push(new Constraint(this, point));
  }

  free(constraint) {   // Removes a constraint, used when the cloth is being cut.
    this.constraints.splice(this.constraints.indexOf(constraint), 1);
  }

  addForce(x, y) {   // Applies a force to the point, affecting its velocity.
    this.vx += x;
    this.vy += y;
  }

  pin(pinx, piny) {  // Pins the point to a specific position. Is being used for fixing the top row of the cloth.
    this.pinX = pinx;
    this.pinY = piny;
  }}


class Constraint {     // Models the link between two points, ensuring they stay within a designated distance or providing spacing.
  constructor(p1, p2) {
    this.p1 = p1;
    this.p2 = p2;
    this.length = spacing;
  }

  resolve() {      //Adjusts points to satisfy the constraint. If the distance exceeds tearDist, the link is broken.

    let dx = this.p1.x - this.p2.x;     // Difference in x positions
    let dy = this.p1.y - this.p2.y;     // Difference in y positions
    let dist = Math.sqrt(dx * dx + dy * dy);    // Distance between points

    // To skip further calculations if the points are already close enough
    if (dist < this.length) return;

    let diff = (this.length - dist) / dist;   // Difference ratio

    // Tear the cloth if the distance is greater than the tear distance
    if (dist > tearDist) this.p1.free(this);

    // Calculate the amount of movement needed to maintain the constraint
    let mul = diff * 0.5 * (1 - this.length / dist);

    let px = dx * mul;
    let py = dy * mul;

    // Move points only if they are not pinned
    !this.p1.pinX && (this.p1.x += px);
    !this.p1.pinY && (this.p1.y += py);
    !this.p2.pinX && (this.p2.x -= px);
    !this.p2.pinY && (this.p2.y -= py);

    return this;
  }

  draw() {   //Draws the line between the two points.

    ctx.moveTo(this.p1.x, this.p1.y);
    ctx.lineTo(this.p2.x, this.p2.y);
  }}


class Cloth {      // Represents the entire cloth made up of points as part of the grid.
  constructor(free) {
    this.points = [];

    let startX = canvas.width / 2 - clothX * spacing / 2;

    for (let y = 0; y <= clothY; y++) {
      for (let x = 0; x <= clothX; x++) {
        let point = new Point(startX + x * spacing, 20 + y * spacing);
        !free && y === 0 && point.pin(point.x, point.y);
        x !== 0 && point.attach(this.points[this.points.length - 1]);
        y !== 0 && point.attach(this.points[x + (y - 1) * (clothX + 1)]);

        this.points.push(point);
      }
    }
  }

  update(delta) {     //Applies physics updates and draws the cloth. 
    let i = 5;

    // Resolve each constraint accuracy times per frame for stability
    while (i--) {      
      this.points.forEach(point => {
        point.resolve();
      });
    }
    //updating and drawing each point.
    ctx.beginPath();
    this.points.forEach(point => {
      point.update(delta * delta).draw();  
    });
    ctx.stroke();   // Finalize the points drawing.
  }}

// Updates the mouse object with the new mouse position relative to the canvas.
function setMouse(e) {
  let rect = canvas.getBoundingClientRect();
  mouse.px = mouse.x;
  mouse.py = mouse.y;
  mouse.x = e.clientX - rect.left;
  mouse.y = e.clientY - rect.top;
}

canvas.onmousedown = e => {
  mouse.button = e.which;
  mouse.down = true;
  setMouse(e);
};

canvas.onmousemove = setMouse;

canvas.onmouseup = () => mouse.down = false;

canvas.oncontextmenu = e => e.preventDefault();

let cloth = new Cloth();

function disableGravity() {
  gravity = 0;
  cloth = new Cloth(true);
}

;(function update(time) {
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  cloth.update(0.016);

  window.requestAnimFrame(update);
})(0);