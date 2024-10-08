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

const N = 300; // number of particles	

const spacing = 1.0; // Spacing of particles

const k = spacing / 1000; // Far pressure weight
//const k = 0.1; // Far pressure weight

const k_near = k*10; // near pressure weight
const gravity = 0.005//
const rest_density = 3;
const r=spacing*1.25;						// Radius of Support
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
    this.sigma = 0.1;
    this.beta = 0.0;
    this.neighbors = [];
  }

    getParticleColorAndSymbol() {	
		// Calculate color components
        const x = 0.1 * this.rho;               // Blue component based on density
        const y = Math.min(1, Math.abs(50 * this.vel.x));
        const z = Math.min(1, Math.abs(50 * this.force.x));  // Red component based on x-velocity


		// let	velocityMagnitude = this.force.len();
		// let	forceMagnitude = this.force.len();
		
		// const y = Math.min(1, Math.abs(50 * velocityMagnitude));
        // const z = Math.min(1, Math.abs(50 * forceMagnitude));  // Red component based on x-velocity


        // Create an RGB color in the range [0, 255]
        const red = Math.floor(255 * (0.3 + x));
        const green = Math.floor(255 * (0.3 + y));
        const blue = Math.floor(255 * (0.9 + z));

        // Determine a symbol based on some velocity magnitude (you can customize this)

		let	velocityMagnitude = this.force.len();
		
        let symbol;

		var fontStyle = "";

        if (velocityMagnitude < 0.003) {
            symbol = 'O';  // Low velocity
			fontStyle = "bold";
        } else if (velocityMagnitude < 0.05) {
            symbol = 'o';  // Medium velocity
			fontStyle = "";
        } else {
            symbol = '*';  // High velocity
			fontStyle = "italic";
        }

        return {
            color: `rgb(${red}, ${green}, ${blue})`, // RGB color string
            symbol: symbol,
			fontStyle: fontStyle
        };
    }
}

function idle() {
	for(let i = 0; i < particles.length; i++) {
		// set the old position to curr position
		particles[i].pos_old = new Vec2(particles[i].pos.x,particles[i].pos.y);

		var forceVelVec = new Vec2(particles[i].vel.x + particles[i].force.x,particles[i].vel.y + particles[i].force.y);

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

		//console.log(particles[i].pos.y);

		if(particles[i].pos.y < 0) {particles[i].force.y -= (particles[i].pos.y - 0) / 8};
		if(particles[i].pos.y > height){particles[i].force.y -= (particles[i].pos.y - height) / 8};


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
		for(let j = 0; j < particles.length; j++) {

			if(j >= i) {
				continue;
			}

			// The vector seperating the two particles
			let rij = new Vec2(particles[j].pos.subVec(particles[i].pos).x,particles[j].pos.subVec(particles[i].pos).y);
			let rij_len2 = rij.lenSquared();

			// If they're within the radius of support ...
			if(rij_len2 < rsq)
			{
				// Get the actual distance from the squared distance.
				let rij_len = Math.sqrt(rij_len2);

				// And calculated the weighted distance values
				let q = 1 - rij_len / r;
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
			let rij = new Vec2( particles[j].pos.x - particles[i].pos.x,particles[j].pos.y - particles[i].pos.y);



			


			// calculate the force from the pressures calculated above
			let dm = (particles[i].press + particles[j].press) * q +
			(particles[i].press_near + particles[j].press_near) * q2;

			// Get the direction of the force
			

			let D = new Vec2((rij.normalize().multScalar(dm)).x,(rij.normalize().multScalar(dm)).y);
			
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
			let rij = new Vec2(particles[n.j].pos.subVec(particles[i].pos).x,particles[n.j].pos.subVec(particles[i].pos).y);
			let l = (rij).len();
			let q = l / r;

			let rijn = new Vec2(rij.divScalar(l));

			

			// Get the projection of the velocities onto the vector between them.
			let u = (particles[n.i].vel.subVec(particles[n.j].vel)).multVec(rijn);
			if(u > 0)
			{
				// Calculate the viscosity impulse between the two particles
				// based on the quadratic function of projected length.
				//let I = new Vec2((1 - q) * (particles[n.j].sigma * u + particles[n.j].beta * u*u) * rijn);
				let I = new Vec2(rijn.multScalar(1 - q).multScalar(particles[n.j].sigma * u + particles[n.j].beta * u*u));

				// Apply the impulses on the two particles
				particles[n.i].vel.subVec(I.multScalar(0.5));
				particles[n.j].vel.addVec(I.multScalar(0.5));
			}

		}
	}	
}









function render() {
    // Create an empty grid
    let grid = Array(height).fill().map(() => Array(width).fill(' '));

    for (const particle of particles) {
        // Convert position to integer for grid placement
        const gridX = Math.floor(particle.pos.x);
        const gridY = Math.floor(particle.pos.y);

        // Make sure particles are within bounds
        if (gridX >= 0 && gridX < width && 
			gridY >= 0 && gridY < height) {
            const { color, symbol, fontStyle } = particle.getParticleColorAndSymbol();

            grid[gridY][gridX] += `<span style="font-style:${fontStyle};color:${color};">${symbol}</span>`;
        }
    }
	
    asciiDisplay.innerHTML = grid.map(row => row.join('')).join('<br>');
}

function init() {
	let w = width / 4;
    for(let y = height - 10; y <= 10000; y += r * 0.5) {
		for(let x = 0; x <= w; x += r * 0.5)
		{
			if(particles.length > N) {
				break;
			}

			let p = new Particle(x,y);
			p.pos_old = p.pos.addVec(new Vec2(0.001,0.001)).multVec(new Vec2(2,3));

			particles.push(p);
		}
	}

}

// Main loop
async function animate() {
    if(firstTime == 0) {
		firstTime = 1;
		init();
	}
	idle();
	render();
    requestAnimationFrame(animate);
}

animate();
