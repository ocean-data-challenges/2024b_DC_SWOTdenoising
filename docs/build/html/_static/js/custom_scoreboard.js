$(function () {
  var includes = $('[data-include]')
  $.each(includes, function () {
    var file = '/Users/sammymetref/Documents/DataChallenges/2024_DC_WOC-ESA/docs/source/_static/' + $(this).data('include') + '.html'
    $(this).load(file)
  })
})