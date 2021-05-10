/* jspsych-instructions.js
 * Josh de Leeuw
 *
 * This plugin displays text (including HTML formatted strings) during the experiment.
 * Use it to show instructions, provide performance feedback, etc...
 *
 * Page numbers can be displayed to help with navigation by setting show_page_number
 * to true.
 *
 * documentation: docs.jspsych.org
 *
 *
 */

jsPsych.plugins.graphs = (function() {

  var plugin = {};

  plugin.info = {
    name: 'graphs_freeSort',
    description: '',
    parameters: {
      pages: {
        type: jsPsych.plugins.parameterType.HTML_STRING,
        pretty_name: 'Pages',
        default: undefined,
        array: true,
        description: 'Each element of the array is the content for a single page.'
      },
      button_label: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button label',
        default:  'Continue',
        description: 'The text that appears on the button to continue to the next trial.'
      }
    }
  }

  plugin.trial = function(display_element, trial) {

    var start_time = performance.now();

    var current_page = 0;

    var view_history = [];

    var start_time = performance.now();

    var last_page_update_time = start_time;

    function btnListener(evt){
    	evt.target.removeEventListener('click', btnListener);
    	if(this.id === "jspsych-instructions-back"){
    		back();
    	}
    	else if(this.id === 'jspsych-instructions-next'){
    		next();
    	}
    }

    function show_current_page() {
      scrollWin();
      var html = trial.pages[current_page];
      
      html += "<p id='accuDiv' width='200'></p>";

      display_element.innerHTML = html;
      
      display_element.innerHTML += '<p><button id="jspsych-free-sort-done-btn" class="jspsych-btn">'+trial.button_label+'</button>';
      
    }

    

    show_current_page();
    Plotly.newPlot('accuDiv', dataPlot, layout, {staticPlot: true});
    
    display_element.querySelector('#jspsych-free-sort-done-btn').addEventListener('click', function(){
      var end_time = performance.now();
      var rt = end_time - start_time;
      var trial_data = {
        "rt": rt,
      };

      // advance to next part
      display_element.innerHTML = '';
      jsPsych.finishTrial(trial_data);
    });

  };
  

  return plugin;
})();
