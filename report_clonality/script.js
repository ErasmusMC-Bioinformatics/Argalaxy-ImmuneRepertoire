function compareAdd(id){
	var img = document.createElement('img');
	img.setAttribute('src', 'HeatmapVD_' + id + '.png');
	var td = document.createElement('td');
	td.setAttribute('id', "comparison_vd_" + id);
	td.appendChild(img)
	$('#comparison_table_vd').append(td);
	
	img = document.createElement('img');
	img.setAttribute('src', 'HeatmapVJ_' + id + '.png');
	td = document.createElement('td');
	td.setAttribute('id', "comparison_vj_" + id);
	td.appendChild(img)
	$('#comparison_table_vj').append(td);
	
	img = document.createElement('img');
	img.setAttribute('src', 'HeatmapDJ_' + id + '.png');
	td = document.createElement('td');
	td.setAttribute('id', "comparison_dj_" + id);
	td.appendChild(img)
	$('#comparison_table_dj').append(td);
	
	$('#compare_checkbox_' + id).attr('onchange', "javascript:compareRemove('" + id + "')");
}


function compareRemove(id){
	$("#comparison_vd_" + id).remove()
	$("#comparison_vj_" + id).remove()
	$("#comparison_dj_" + id).remove()
	$("#compare_checkbox_" + id).attr('onchange', "javascript:compareAdd('" + id + "')");
}

$( document ).ready(function () {
	$('#junction_table').tablesorter();
})
