class Neighbor {
constructor(i = 0, j = 0, q = 0, q2 = 0) {
	this.i = i;
	this.j = j;

	this.q = q;
	this.q2 = q2;
}
};

class Vec2 {
  constructor(x = 0, y = 0) {
    this.x = x;
    this.y = y;
  }

  addVec(b) {
    return new Vec2(this.x + b.x, this.y + b.y);
  }

  addScalar(scalar) {
    return new Vec2(this.x + scalar, this.y + scalar);
  }

  subVec(b) {
    return new Vec2(this.x - b.x, this.y - b.y);
  }

  subScalar(scalar) {
    return new Vec2(this.x - scalar, this.y - scalar);
  }

  multVec(b) {
    return new Vec2(this.x * b.x, this.y * b.y);
  }

  multScalar(scalar) {
    return new Vec2(this.x * scalar, this.y * scalar);
  }

  divScalar(scalar) {
    return new Vec2(this.x / scalar, this.y / scalar);
  }

  lenSquared() {
    return this.x * this.x + this.y * this.y;
  }

  len() {
    return Math.sqrt(this.lenSquared());
  }

  normalize() {
    return this.divScalar(this.len());
  }
}

function clamp(number, min, max) {
    return Math.max(min, Math.min(number, max));
  }

const N = 500; // number of particles	

const spacing = 1.0; // Spacing of particles

const k = spacing / 1000; // Far pressure weight
//const k = 0.1; // Far pressure weight

const k_near = k*2; // near pressure weight
const gravity = 0.005 //
const rest_density = 3;
const r=spacing*1.25;						// Radius of Support
//const r=spacing;						// Radius of Support
const rsq=r*r;									// ... squared for performance stuff

const width = 100;
const height = 30;
const particles = [];


const asciiDisplay = document.getElementById('ascii-render');

var firstTime = 0;

class Particle {
  constructor(x, y) {
    this.pos = new Vec2(x, y);
    this.pos_old = new Vec2(x, y);
    this.vel = new Vec2(0, 0);
    this.force = new Vec2(0, 0);
    this.mass = 1.0;
    this.rho = 0;
    this.rho_near = 0;
    this.press = 0;
    this.press_near = 0;
	
    // if a highly viscous behavior is desired, σcan be increased. 
	// For less viscous ﬂuids, only β should be set to a non-zero value
	this.sigma = 0.1;
    this.beta = 0.0;

    this.neighbors = [];
  }

    getParticleColorAndSymbol() {	
		// Calculate color components
        const x = 0.1 * this.rho;               // Blue component based on density

		//console.log("rho",this.vel.x);

        const y = Math.min(1, Math.abs(50 * this.vel.x));
        const z = Math.min(1, Math.abs(50 * this.force.x));  // Red component based on x-velocity


		// let	velocityMagnitude = this.force.len();
		let	forceMagnitude = this.force.len();
		
		// const y = Math.min(1, Math.abs(50 * velocityMagnitude));
        // const z = Math.min(1, Math.abs(50 * forceMagnitude));  // Red component based on x-velocity


		//let r = 255 * proximityToPoint(this.vel.x, 0.3);
		//let g = 255 * proximityToPoint(this.vel.y,0.6);
		//let b = 255 * proximityToPoint(forceMagnitude,0.2);

		let scale = proximityToPoint(this.vel.y,0.6) *100;

		let color1 = {r:189,g:16,b:63};
		let color2 = {r:191,g:21,b:92};

		

		const colors = colorInterpolate(color1,color2,scale)
		//console.log(colors);

		//let r = 255;
		//let g = 255 * proximityToPoint(this.vel.y,0.6);
		//let b = 100;
		//let b = 255 * proximityToPoint(this.vel.x,0.6);

		//console.log(red,green,blue);
		//const colors = scaleColorsTo255(r,g,b);

        // Create an RGB color in the range [0, 255]
        //const red = Math.floor(255 * (0.3 + x));
        //const green = Math.floor(255 * (0.3 + y));
       // const blue = Math.floor(255 * (0.9 + z));

        // Determine a symbol based on some velocity magnitude (you can customize this)

		let	velocityMagnitude = this.force.len();
		
        let symbol;

		var fontStyle = "";

		console.log(scale);

		let colorMag = (colors.r + colors.g + colors.b) / 765.0;

        if (velocityMagnitude < 0.003) {
            symbol = 'O';  // Low velocity
			fontStyle = "bold";
        } else if (velocityMagnitude < 0.006) {
            symbol = 'o';  // Medium velocity
			fontStyle = "";
        } else {
            symbol = '0';  // High velocity
			fontStyle = "italic";
        }

        // if (velocityMagnitude < 0.003) {
        //     symbol = 'O';  // Low velocity
		// 	fontStyle = "bold";
        // } else if (velocityMagnitude < 0.04) {
        //     symbol = 'X';  // Medium velocity
		// 	fontStyle = "";
        // } else {
        //     symbol = '*';  // High velocity
		// 	fontStyle = "italic";
        // }



        return {
            color: `rgb(${colors.r}, ${colors.g}, ${colors.b})`, // RGB color string
            symbol: symbol,
			fontStyle: fontStyle
        };
    }
}

function colorInterpolate(colorA, colorB, intval) {
	return {
	r : Math.round(colorA.r * (1 - intval) + colorB.r * intval),
	g : Math.round(colorA.g * (1 - intval) + colorB.g * intval),
	b : Math.round(colorA.b * (1 - intval) + colorB.b * intval)
	};
  }

function scaleColorsTo255(r, g, b) {
    // Find the maximum value among x, y, and z
    const maxVal = Math.max(r, g, b);

    // If the maximum value is 0, scaling isn't possible (return the original values)
    if (maxVal === 0) {
        return { r, g, b };
    }

    // Calculate the scaling factor
    const scaleFactor = 255 / maxVal;

    // Scale each variable
    const scaledR = Math.min(Math.round(r * scaleFactor), 255);
    const scaledG = Math.min(Math.round(g * scaleFactor), 255);
    const scaledB = Math.min(Math.round(b * scaleFactor), 255);

    return { r: scaledR, g: scaledG, b: scaledB };
}

function proximityToPoint(x, target) {
    //const point = 0.2;  // The point of interest
    //const maxDistance = 0.2;  // Maximum distance we care about (beyond this, output will be close to 0.01)

    // Calculate the distance from x to the target
    const distance = Math.abs(x - target);

    // Normalize the distance between 0 and 1
    const normalizedDistance = Math.min(distance / target, 1);

    // Invert the normalized distance to get a higher value when closer
    const invertedDistance = 1 - normalizedDistance;

    // Scale the inverted distance to the range [0.01, 1]
    const result = 0.01 + (invertedDistance * (1 - 0.01));

    return result;
}


function idle() {
	for(let i = 0; i < particles.length; i++) {
		// set the old position to curr position
		particles[i].pos_old = new Vec2(particles[i].pos.x,particles[i].pos.y);

		var forceVelVec = particles[i].vel.addVec(particles[i].force);

		// add velocity and force to the current position
		particles[i].pos = particles[i].pos.addVec(forceVelVec);

		// Restart the forces with gravity only. We'll add the rest later.
		particles[i].force = new Vec2(0,gravity);

		// Calculate the velocity for later.
		particles[i].vel = particles[i].pos.subVec(particles[i].pos_old);
		// If the velocity is really high, we're going to cheat and cap it.
		// This will not damp all motion. It's not physically-based at all. Just
		// a little bit of a hack.
		let max_vel = 2.0;

		let vel_mag = particles[i].vel.lenSquared();
		// If the velocity is greater than the max velocity, then cut it in half.
		if(vel_mag > max_vel * max_vel) {

			particles[i].vel = particles[i].vel.multScalar(0.5);
		}
		
		// If the particle is outside the bounds of the world, then
		// Make a little spring force to push it back in.
		if(particles[i].pos.x < 0) { particles[i].force.x -= (particles[i].pos.x - 0) / 8};
		if(particles[i].pos.x >  width) {particles[i].force.x -= (particles[i].pos.x - width) / 8};

		if(particles[i].pos.y < 0) {particles[i].force.y -= (particles[i].pos.y - 0) / 8};
		
		let temp = particles[i].force.y;
		if(particles[i].pos.y > height){particles[i].force.y -= (particles[i].pos.y - height) / 16};
		//console.log(temp,particles[i].force.y);

		// Reset the nessecary items.
		particles[i].rho = 0; 
		particles[i].rho_near = 0;
		particles[i].neighbors = [];
	}

	// TODO: call calc density, pressure, etc

	CalcDensity();
	CalcPressureAndForce();
	CalcViscosity();

}

// Calculate the density by basically making a weighted sum
// of the distances of neighboring particles within the radius of support (r)
function CalcDensity() {
	for(let i = 0; i < particles.length; i++) { // for each particle
		particles[i].rho = 0; 
		particles[i].rho_near = 0;

		// We will sum up the 'near' and 'far' densities.
		let d = 0;
		let dn = 0;

		// only look at each pair of particles once. 
		// dont calc interaction for a particle with itself!
		for(let j = 0; j < i; j++) {

			// The vector seperating the two particles
			let rij = particles[j].pos.subVec(particles[i].pos);

			let rij_len2 = rij.lenSquared();

			// If they're within the radius of support ...
			if(rij_len2 < rsq)
			{
				// Get the actual distance from the squared distance.
				let rij_len = Math.sqrt(rij_len2);

				// And calculated the weighted distance values
				let q = (1 - rij_len) / r;
				let q2 = q*q;
				let q3 = q2*q;

				d += q2;
				dn += q3;

				// Accumulate on the neighbor
				particles[j].rho += q2;
				particles[j].rho_near += q3;

				// Set up the neighbor list for faster access later.
				let n = new Neighbor(i,j,q,q2);        
				particles[i].neighbors.push(n);
			}
		}
		particles[i].rho += d;
		particles[i].rho_near += dn;
	}
}

function CalcPressureAndForce() {
	// Make the simple pressure calculation from the equation of state.
	
	for(let i = 0; i < particles.length; i++)
	{
		particles[i].press = k * (particles[i].rho - rest_density);
		particles[i].press_near = k_near * particles[i].rho_near;
	}

	// Force particles in or out from their neighbors based on their difference from the rest density

	for(let i = 0; i < particles.length; i++)
	{
		var dX = new Vec2(0,0);

		// For each of the neighbors
		let ncount = particles[i].neighbors.length;
		
		for(let ni = 0; ni < ncount; ni++)
		{
			//let n = particles[i].neighbors[ni]; // TODO?
			let j = particles[i].neighbors[ni].j;     
			
			let q = particles[i].neighbors[ni].q;
			let q2 = particles[i].neighbors[ni].q2;


			// The vector from particle i to particle j
			let rij = particles[j].pos.subVec(particles[i].pos); 

			// calculate the force from the pressures calculated above
			let dm = (particles[i].press + particles[j].press) * q +
			(particles[i].press_near + particles[j].press_near) * q2;

			// Get the direction of the force
			

			//let D = new Vec2((rij.normalize().multScalar(dm)).x,(rij.normalize().multScalar(dm)).y);
			let D = rij.normalize().multScalar(dm);
			
			//			dX = new Vec2(dX.addVec(D).x,dX.addVec(D).y);
			dX = dX.addVec(D);

			particles[j].force = particles[j].force.addVec(D);

			
		}
		particles[i].force.subVec(dX);
		
	}
}

function CalcViscosity() {
	for(let i=0; i < particles.length; i++)
	{
		// For each of that particles neighbors
		for(let ni=0; ni < particles[i].neighbors.length; ni++)
		{
			let n = particles[i].neighbors[ni]; // TODO?

			// The vector from particle i to particle j
			let rij = particles[n.j].pos.subVec(particles[i].pos); 
			
			let l = (rij).len();
			let q = l / r;

			let rijn = rij.divScalar(l);

			

			// Get the projection of the velocities onto the vector between them.
			let u = (particles[n.i].vel.subVec(particles[n.j].vel)).multVec(rijn);
			if(u > 0)
			{
				// Calculate the viscosity impulse between the two particles
				// based on the quadratic function of projected length.

				//let temp = (1 - q) * (particles[n.j].sigma * u + particles[n.j].beta * u*u);
				//let I = new Vec2(rijn.x * temp,rijn.y*temp);
				
				let I = rijn.multScalar( (1 - q)) * (particles[n.j].sigma * u + particles[n.j].beta * u*u);

				// Apply the impulses on the two particles
				particles[n.i].vel.subVec(I.multScalar(0.5));
				particles[n.j].vel.addVec(I.multScalar(0.5));
			}

		}
	}	
}


// Persistent grid to avoid recreating it each frame
let grid = Array(height).fill().map(() => Array(width).fill(' '));

function render() {
    // Clear the grid in a more optimized way
    for (let y = 0; y < height; y++) {
        grid[y].fill(' '); // Reset each row
    }

    let html = ''; // Batch the html string

    for (const particle of particles) {
        var gridX = Math.floor(particle.pos.x);
        var gridY = Math.floor(particle.pos.y);

		gridX = clamp(gridX,0,width);
		gridY = clamp(gridY,0,height);

        if (gridX >= 0 && gridX < width && gridY >= 0 && gridY < height) {
            const { color, symbol, fontStyle } = particle.getParticleColorAndSymbol();
            // Directly update the grid's cell with the HTML string for this particle
            grid[gridY][gridX] = `<span style="font-style:${fontStyle};color:${color};">${symbol}</span>`;
        } else {
			//console.log(gridX,gridY);
		}
    }

    // Build up the HTML string outside the loop
    html = grid.map(row => row.join('')).join('<br>');

    // Only update the DOM if necessary
    if (asciiDisplay.innerHTML !== html) {
        asciiDisplay.innerHTML = html;
    }
}



// function render() {
//     // Create an empty grid
//     let grid = Array(height).fill().map(() => Array(width).fill(' '));

//     for (const particle of particles) {
//         // Convert position to integer for grid placement
//         const gridX = Math.floor(particle.pos.x);
//         const gridY = Math.floor(particle.pos.y);

//         // Make sure particles are within bounds
//         if (gridX >= 0 && gridX < width && 
// 			gridY >= 0 && gridY < height) {
//             const { color, symbol, fontStyle } = particle.getParticleColorAndSymbol();

//             grid[gridY][gridX] = `<span style="font-style:${fontStyle};color:${color};">${symbol}</span>`;
//             // grid[gridY][gridX] += `<span style="font-style:${fontStyle};color:${color};">${symbol}</span>`;
//         }
//     }
	
//     asciiDisplay.innerHTML = grid.map(row => row.join('')).join('<br>');
// }

function init() {
	let w = width / 4;
    for(let y = height - 10; y <= 10000; y += r * 0.5) {
		for(let x = 0; x <= w; x += r * 0.5)
		{
			if(particles.length > N) {
				break;
			}

			let p = new Particle(x,y);
			p.pos_old = p.pos.addVec(new Vec2(0,0)).multVec(new Vec2(2,3));

			particles.push(p);
		}
	}

}





var stop = false;
var frameCount = 0;
var fps, fpsInterval, startTime, now, then, elapsed;


// initialize the timer variables and start the animation

function startAnimating(fps) {
    fpsInterval = 1000 / fps;
    then = Date.now();
    startTime = then;
    animate();
}



// Main loop
// async function animate() {
//     if(firstTime == 0) {
// 		firstTime = 1;
// 		init();
// 	}
// 	idle();
// 	render();
//     requestAnimationFrame(animate);
// }
// animate();

startAnimating(30);

// the animation loop calculates time elapsed since the last loop
// and only draws if your specified fps interval is achieved
function animate() {
    // request another frame

	if(firstTime == 0) {
		firstTime = 1;
		init();
	}

	idle(); // always idle?
	//idle(); // always idle?

    requestAnimationFrame(animate);

    // calc elapsed time since last loop

    now = Date.now();
    elapsed = now - then;

    // if enough time has elapsed, draw the next frame

    if (elapsed > fpsInterval) {

        // Get ready for next frame by setting then=now, but also adjust for your
        // specified fpsInterval not being a multiple of RAF's interval (16.7ms)
        then = now - (elapsed % fpsInterval);

        // drawing code here
		render();
    }
}