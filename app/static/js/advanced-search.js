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