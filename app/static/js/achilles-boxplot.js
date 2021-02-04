function plotAchillesBoxplot(bicluster, perRow, boxplotSizePerDatum) {
	d3.select("#achilles-boxplot svg").remove()
	
	const main = d3.select("#achilles-boxplot")
		.append("svg")
		.attr("id", "boxplot")
		.style("background-color", "#FFF")
	
	const boxWidth = perRow * boxplotSizePerDatum
	const boxHeight = 400
	const heightMargin = 60

	const tooltip = d3.select("#achilles-tooltip")
	
	d3.json(`/achilles-bicluster/${bicluster}`, (json) => {
		json = json.sort((a, b) => a.name.localeCompare(b.name))
		main.attr("width", boxWidth + 50)
			.attr("height", 10 + (boxHeight + heightMargin) * (Math.ceil(json.length / perRow)))
		
		for(let i = 0; i < json.length; i += perRow) {
			const svg = main.append("g")
				.attr("transform", `translate(50, ${10 + (boxHeight + heightMargin) * (i / perRow)})`)
			
			const slice = json.slice(i, i + perRow)

			// add separation line
			main.append("line")
				.attr("x1", 0)
				.attr("x2", boxWidth + 50)
				.attr("y1", 10 + (boxHeight + heightMargin) * (i / perRow) - heightMargin / 2 + 10)
				.attr("y2", 10 + (boxHeight + heightMargin) * (i / perRow) - heightMargin / 2 + 10)
				.attr("stroke", "#333")
				.attr("stroke-width", "2px")
			
			// x-axis
			const x = d3.scaleBand()
				.range([0, boxWidth])
				.domain(slice.map(a => a.name))
				.paddingInner(1)
				.paddingOuter(0.5)
		
			svg.append("g")
				.attr("class", `xAxis${i}`)
				.attr("transform", "translate(0," + boxHeight + ")")
				.call(d3.axisBottom(x))
			
			// x-axis styling
			d3.select(`.xAxis${i} path`)
				.style("stroke", "transparent")
			
			d3.selectAll(`.xAxis${i} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.xAxis${i} .tick text`)
				.style("font-size", "11pt")
			
			// y-axis
			const minY = -2.1, maxY = 2.1
			const y = d3.scaleLinear()
				.domain([minY, maxY])
				.range([boxHeight, 0])

			svg.append("g")
				.attr("class", `yAxis${i}`)
				.call(d3.axisLeft(y).ticks(8).tickSize(-boxWidth))

			// y-axis styling
			d3.select(`.yAxis${i} path`)
				.style("stroke-width", "0px")
			
			d3.selectAll(`.yAxis${i} line`)
				.style("stroke", "#AAA")
			
			d3.selectAll(`.yAxis${i} .tick`)
				.style("stroke-width", "2px")
			
			d3.selectAll(`.yAxis${i} .tick text`)
				.style("font-size", "10pt")
			
			const tooltipify = (o) => {
				return o.on("mouseover", d => {
					tooltip.transition()
						.duration(200)
						.style("opacity", 0.9)
					
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
					
					tooltip.html(`<b>${d.name}</b><br />Maximum: ${d.high.toFixed(4)}<br />Upper Quartile: ${d.q3.toFixed(4)}<br />Median: ${d.median.toFixed(4)}<br />Loewr Quartile: ${d.q1.toFixed(4)}<br />Minimum: ${d.low.toFixed(4)}`)
				})
				.on("mousemove", d => {
					tooltip.style("left", `${d3.event.clientX}px`)
						.style("top", `${d3.event.clientY}px`)
				})
				.on("mouseout", d => {
					tooltip.transition()
						.duration(200)
						.style("opacity", 0)
					
					tooltip.style("left", `-1000px`)
						.style("top", `-1000px`)
				})
			}

			const length = 50
			// line through -1
			svg.append("line")
				.attr("x1", 0)
				.attr("x2", boxWidth + 50)
				.attr("y1", Math.floor(y(-1)))
				.attr("y2", Math.floor(y(-1)))
				.attr("stroke", "red")
				.attr("stroke-width", "2px")

			// enter data
			tooltipify(
				svg.selectAll("boxplot")
					.data(slice)
					.enter()
					.append("line")
					.attr("x1", d => Math.floor(x(d.name)))
					.attr("x2", d => Math.floor(x(d.name)))
					.attr("y1", d => Math.floor(y(d.low)))
					.attr("y2", d => Math.floor(y(d.high)))
					.attr("stroke", "black")
					.attr("stroke-width", "2px")
			)

			// max horizontal line
			tooltipify(
				svg.selectAll("boxplot")
					.data(slice)
					.enter()
					.append("line")
					.attr("x1", d => Math.floor(x(d.name) - length))
					.attr("x2", d => Math.floor(x(d.name) + length))
					.attr("y1", d => Math.floor(y(d.high)))
					.attr("y2", d => Math.floor(y(d.high)))
					.attr("stroke", "black")
					.attr("stroke-width", "2px")
			)
			
			// min horizontal line
			tooltipify(
				svg.selectAll("boxplot")
					.data(slice)
					.enter()
					.append("line")
					.attr("x1", d => Math.floor(x(d.name) - length))
					.attr("x2", d => Math.floor(x(d.name) + length))
					.attr("y1", d => Math.floor(y(d.low)))
					.attr("y2", d => Math.floor(y(d.low)))
					.attr("stroke", "black")
					.attr("stroke-width", "2px")
			)
			
			// drawing the boxes
			tooltipify(
				svg.selectAll("boxplot")
					.data(slice)
					.enter()
					.append("rect")
					.attr("x", d => Math.floor(x(d.name) - length))
					.attr("y", d => Math.floor(y(d.q3)))
					.attr("height", d => Math.floor(y(d.q1) - y(d.q3)))
					.attr("width", length * 2)
					.attr("stroke", "black")
					.attr("stroke-width", "2px")
					.style("fill", d => "#E6F4FF")
			)
			
			// median line
			tooltipify(
				svg.selectAll("boxplot")
					.data(slice)
					.enter()
					.append("line")
					.attr("x1", d => Math.floor(x(d.name) - length))
					.attr("x2", d => Math.floor(x(d.name) + length))
					.attr("y1", d => Math.floor(y(d.median)))
					.attr("y2", d => Math.floor(y(d.median)))
					.attr("stroke", "black")
					.attr("stroke-width", "2px")
			)
		}
	})
}