function addFilterSelect(id, options) {
	let div = document.createElement("div")
	let trash = document.createElement("span")
	trash.className = "glyphicon glyphicon-trash"
	trash.addEventListener("click", (event) => {
		div.remove()
	})

	let select = document.createElement("select")
	select.className = "form-control"
	select.value = options[0]
	for(let option of options) {
		let optionNode = document.createElement("option")
		optionNode.value = option
		optionNode.innerHTML = option
		select.append(optionNode)
	}

	div.append(trash, select)
	document.getElementById(id).append(div)
}

function createFilter(buttonId, selectId, containerId, options) {
	document.getElementById(buttonId).addEventListener("click", (event) => {
		addFilterSelect(containerId, options)
	})

	document.getElementById(selectId).value = options[0]
	for(let option of options) {
		let optionNode = document.createElement("option")
		optionNode.value = option
		optionNode.innerHTML = option
		document.getElementById(selectId).append(optionNode)
	}
}

const idToName = {
	1: "Self sufficiency in growth signals",
	2: "Reprogramming energy metabolism",
	3: "Evading apoptosis",
	4: "Genome instability and mutation",
	5: "Sustained angiogenesis",
	6: "Tissue invasion and metastasis",
	7: "Tumor promoting inflammation",
	8: "Limitless replicative potential",
	9: "Evading immune detection",
	10: "Insensitivity to antigrowth signals",
}

const contexts = new Set()
const activeContexts = new Set()
const indexToName = {}
function addHallmarkImage(index, grayscale, context) {
	let width = 300, height = 274
	if(!context) {
		let container = document.getElementById("hallmark-container")
		let canvas = document.createElement("canvas")
		canvas.width = width
		canvas.height = height
		container.appendChild(canvas)
		context = canvas.getContext("2d")
	}

	if(grayscale) {
		context.filter = "grayscale(100%) opacity(30%)"	
	}
	else {
		context.filter = "none"
	}

	context.clearRect(0, 0, width, height)
	context.drawImage(document.getElementById(`hallmark-${index}`), 0, 0, width, height)
	context.imageSource = index

	if(!contexts.has(context)) {
		contexts.add(context)
	}
}

const selectedHallmarkNames = []
function onHallmarkClick(event) {
	let foundContext
	for(let context of contexts) {
		let rgba = context.getImageData(event.offsetX, event.offsetY, 1, 1).data
		if(rgba.reduce((sum, curr) => sum + curr) != 0) {
			foundContext = context
		}
	}
	
	if(foundContext) {
		addHallmarkImage(foundContext.imageSource, activeContexts.has(foundContext), foundContext)

		if(activeContexts.has(foundContext)) { // we've deselected the hallmark
			activeContexts.delete(foundContext)
			selectedHallmarkNames.splice(selectedHallmarkNames.indexOf(foundContext.imageSource), 1)
		}
		else { // we've selected the hallmark
			activeContexts.add(foundContext)
			selectedHallmarkNames.push(foundContext.imageSource)
		}

		document.getElementById("hallmark").value = btoa(JSON.stringify(selectedHallmarkNames))
	}
}

document.getElementById("hallmark-container").addEventListener("click", event => onHallmarkClick(event))

let loadedImages = 0
function loadHallmarkImage() {
	loadedImages++
	if(loadedImages == 10) {
		for(let i = 1; i <= 10; i++) {
			addHallmarkImage(i, true)
		}
	}
}

for(let i = 1; i <= 10; i++) {
	$(`#hallmark-${i}`).load(() => loadHallmarkImage())
}