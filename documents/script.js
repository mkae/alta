function Compare(id, width, height) {
	if (!(this instanceof Compare)) {
		return new Compare(id, width, height);
	}

	var canvas = document.createElement('canvas'),
		 container = document.getElementById(id),
		 divide = 0.5;

	this.id = id;
	this.ctx = canvas.getContext('2d');
	this.images = [];
	this.texts  = [];

	// Event handlers
	canvas.addEventListener('mousemove', handler, false);
	canvas.addEventListener('mousedown', handler, false);
	canvas.addEventListener('mouseup', handler, false);

	var self = this;

	function handler(ev) {
		
		ev._x = ev.pageX - this.offsetLeft;
		ev._y = ev.pageY - this.offsetBottom;

		var eventHandler = self[ev.type];
		if (typeof eventHandler == 'function') {
			eventHandler.call(self, ev);
		}
	}

	// Draw canvas into its container
	canvas.setAttribute('width', width);
	canvas.setAttribute('height', height);
	container.appendChild(canvas);

	var rect = container.getBoundingClientRect();
	this.posX = rect.left;
	this.width = rect.right - rect.left;

	Object.defineProperty(this, 'ready', {
		get: function () {
			return this.images.length >= 2;
		}
	});

	Object.defineProperty(this, 'width', {
		get: function () {
			return width;
		}
	});

	Object.defineProperty(this, 'height', {
		get: function () {
			return height;
		}
	});

	Object.defineProperty(this, 'divide', {
		get: function () {
			return divide;
		},
		set: function (value) {
			if (value > 1) {
				value = (value / 100);
			}

			divide = value;
			this.draw();
		}
	});
}

Compare.prototype = {
	add: function (src, text) {
		var img = createImage(src, onload.bind(this));

		function onload(event) {
			this.images.push(img);
			this.texts.push(text);

			if (this.ready) {
				this.draw();
			}
		}
	},

	draw: function () {
		if (!this.ready) {
			return;
		}

		var lastIndex = 1;/*this.images.length - 1,*/
			 before = this.images[lastIndex - 1],
			 after = this.images[lastIndex];

		this.drawImages(this.ctx, before, after);
		this.drawHandle(this.ctx);

		var text = this.texts[lastIndex - 1];
		var text_width = this.ctx.measureText(text).width;

		var split = this.divide * this.width;
		this.ctx.fillText(this.texts[lastIndex], split+10, this.height-10);
		this.ctx.fillText(text, split-10-text_width, this.height-10);
	},
	
	update: function () {
		var container = document.getElementById(this.id);
		var rect = container.getBoundingClientRect();
		this.posX = rect.left;
		this.width = rect.right - rect.left;

		console.log(this.id + ' posX = ' + this.posX);
	},

	drawImages: function (ctx, before, after) {
		var split = this.divide * this.width;

		ctx.drawImage(after, 0, 0);
		ctx.drawImage(before, 0, 0, split, this.height, 0, 0, split, this.height);
	},

	drawHandle: function (ctx) {
		var split = this.divide * this.width;

		ctx.fillStyle = "rgb(220, 50, 50)";
		ctx.fillRect(split - 1, 0, 2, this.height);
	},

	mousedown: function (event) {
		var divide = (event._x - this.posX) / this.width;
		//var divide = event._x / this.width;
		this.divide = divide;

		this.dragstart = true;
	},

	mousemove: function (event) {
		if (this.dragstart === true) {
			var divide = (event._x - this.posX) / this.width;
			//var divide = event._x / this.width;
			this.divide = divide;
		}
	},

	mouseup: function (event) {
		var divide = (event._x - this.posX) / this.width;
		//var divide = event._x / this.width;
		this.divide = divide;

		this.dragstart = false;
	}
};




function createImage(src, onload) {
	var img = document.createElement('img');
	img.src = src;

	if (typeof onload == 'function') {
		img.addEventListener('load', onload);
	}

	return img;
}


//var compare = Compare('compare-1', 150, 149);
//compare.add('sample_result_03.png', 'MERL');
//compare.add('sample_result_03.png', 'ALTA, Blinn BRDF');
