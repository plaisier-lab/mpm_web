(function() {
	const sideNav = document.createElement("div")
	sideNav.id = "side-nav"
	const elements = $(".bs-callout-default")
	const names = elements.map((index, element) => element.getAttribute("id")).get()
	const sideNavLinks = {}
	for(const name of names) {
		const element = document.createElement("a")
		element.innerHTML = name
		element.href = `#${name}`
		sideNav.appendChild(element)

		sideNavLinks[name] = element
	}

	document.body.appendChild(sideNav)

	let lastOnScreen
	const checkPosition = () => {
		// check element positions
		const positions = elements.map(
			(index, element) => ({
				name: element.getAttribute("id"),
				position: element.getBoundingClientRect().bottom - 70,
			})
		).get()

		// get the element that is the most on screen
		for(const { name, position } of positions) {
			if(position > 0) {
				sideNavLinks[name].classList.add("active")

				if(lastOnScreen && lastOnScreen !== sideNavLinks[name]) {
					lastOnScreen.classList.remove("active")
				}

				lastOnScreen = sideNavLinks[name]
				break
			}
		}

		window.requestAnimationFrame(checkPosition)
	}
	checkPosition()
})()