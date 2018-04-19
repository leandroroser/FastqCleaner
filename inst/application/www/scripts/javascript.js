
$( document ).ready(function() {


$('html').attr('lang', 'en');


$("#navbar").append("<a href = '../help/index.html' target='_blank'><img id='help' src = 'images/help.svg'></a>");


       var screenWidth = $(window).width();
       $('#navbar').css('width', screenWidth + 'px');
       $('#navbar').css('min-width', screenWidth + 'px');


var acc = document.getElementsByClassName("accordion");
var i;

for (i = 0; i < acc.length; i++) {
    acc[i].onclick = function(){
        theOffset = $(this).offset();
        $('body,html').animate({ 
            scrollTop: theOffset.top - 125
        }, 300);
        this.classList.toggle("active");
        this.nextElementSibling.classList.toggle("show");
        
  };
}


 $("#advancedMenu").hide();
 
$("#advancedButton").click(function(){
    $("#advancedMenu").toggle();
});



 $("#advancedAdapterMenu").hide();
 
$("#advancedAdapterButton").click(function(){
    $("#advancedAdapterMenu").toggle();
});




$('.navbar').addClass('nav-down');
var didScroll;
var lastScrollTop = 0;
var delta = 5;
var navbarHeight = $('.navbar').outerHeight();

$(window).scroll(function(event){
    didScroll = true;
});

setInterval(function() {
    if (didScroll) {
        hasScrolled();
        didScroll = false;
    }
}, 250);

function hasScrolled() {
    var st = $(this).scrollTop();
    
    if(Math.abs(lastScrollTop - st) <= delta)
        return;
    
    if (st > lastScrollTop && st > navbarHeight){
        $('.navbar').removeClass('nav-down').addClass('nav-up');
    } else {
    if(st + $(window).height() < $(document).height() && st < navbarHeight/2 ) {
         $('.navbar').removeClass('nav-up').addClass('nav-down');
        }
    }
    
    lastScrollTop = st;
}


$('#distPlot').load(function() {
  $('cssload-internal').fadeOut();
});
    
});


