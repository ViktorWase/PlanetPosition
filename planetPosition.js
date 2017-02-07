/**
 * This code uses data provided by JPL to calculate the approximate positions
 * of all of the major planets in out solar system.
 */

/** @namespace */
var planetPosition = {

  /**
   *  Solves x=y-c*sin(y) numerically
   *
   *  @param {number} x - the input variable.
   *
   *  @param {number} c - the elliptical constant (in radians).
   *
   *  @return {number} the solution to Kepler's equation.
   */
  keplersEq: function(x, c){
      var tol = 1.0e-6;
      var diff = 100000;
      var y = x + Math.sin(x);
      var d_x, d_y;
      while(Math.abs(diff)>tol){
          d_x = x - (y-c*Math.sin(y));
          d_y = d_x/(1-c*Math.cos(y));
          y += d_y;
          diff = d_y; //This is the diff of the previous iteration, but still.
      }
      return y;
  },


  /**
   *  Converts from unix time (in milliseconds) to modified Julian Date.
   *
   *  @param {number} unixMillisecs - milliseconds since Unix epoch.
   *
   *  @return {number} The modified Julian Date.
   */
  getJulianFromUnix: function(unixMillisecs){
      return unixMillisecs / 86400000 + 2440587.5;
  },


  /**
   *  Calculates the position of a planet. The position is measured in AU (the mean
   *  distance from the sun to the earth) and the time in years.
   *  The position is given in an ecliptic coordinate system, and t=2000 means
   *  the turn of the latest millenium.
   *
   *  @param {number|string} planet - the name or index of the planet. The index should be 0-indexed.
   *
   *  @param {number} t - time measured in years since the birth of Christ. (or milliseconds since
   *    1 Jan 1970 0:00 if unix_time is used)
   *
   *  @param {dictionary} options - the options should be sent in a dictionary. For example:
   *  {icrf:true, unix_time:true} if one wishes to calculate using unix time and
   *  the icrf coordinate systems.
   *  - circular:true uses an approximation of the elliptical orbit using a circle. This
   *  is slightly faster.
   *  - unix_time:true means that the input time is given in illiseconds since 1 Jan 1970 0:00.
   *  - icrf:true changes the coordinat system to an equatoral coordinate system.
   *
   *  @return {Array} - An array of length 3, with the position of the planet.
   */
  getPos: function(planet, t, options){

      //Set default values.
      var unix_time, icrf, circular;
      if (options === undefined) {
        unix_time = false;
        circular = false;
        icrf = false;
      }
      else{
        if (options[unix_time] === undefined){
          unix_time = false;
        }
        else{
          unix_time = true;
        }

        if (options[icrf] === undefined){
          icrf = false;
        }
        else{
          icrf = true;
        }

        if (options[circular] === undefined){
          circular = false;
        }
        else{
          circular = true;
        }
      }

      //The Keplerian parameters and their derivatives.
      var parameters = [[0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593, 0.00000037,
          0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081], [0.72333566, 0.00677672,
          3.39467605, 181.97909950, 131.60246718, 76.67984255, 0.00000390, -0.00004107, -0.00078890,
          58517.81538729, 0.00268329, -0.27769418], [1.00000261, 0.01671123, -0.00001531, 100.46457166,
          102.93768193, 0.0, 0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0],
          [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891, 0.00001847,
          0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343], [5.20288700, 0.04838624,
          1.30439695, 34.39644051, 14.72847983, 100.47390909, -0.00011607, -0.00013253, -0.00183714,
          3034.74612775, 0.21252668, 0.20469106],[9.53667594, 0.05386179, 2.48599187, 49.95424423,
          92.59887831, 113.66242448, -0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216,
          -0.28867794], [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503,
          -0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589], [30.06992276,
          0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574, 0.00026291, 0.00005105,
          0.00035372, 218.45945325, -0.32241464, -0.00508664]];
      var planet_names = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune'];
      var planet_name;
      var planet_idx=2;
      if (typeof(planet) == 'number'){
        planet_idx = Math.round(planet);
        if(Math.abs(planet_idx-planet)>0.0001){
          throw Error("The planet index cannot be a fraction!");
        }
        if (planet_idx<0 || planet_idx>=8){
          throw Error("The planet is out of range! It should be 0-indexed.");
        }
      }
      else if( typeof(planet)=='string'){
        var has_found_it = false;
        for(var i=0; i<planet_names.length; i++){
          if(planet == planet_names[i]){
            has_found_it = true;
            planet_idx = i;
            break;
          }
        }
        if(has_found_it === false){
          throw Error("The planet argument has to be the name of an actual planet!");
        }
      }
      else{
        throw Error("The planet argument has to be an integer or a string!");
      }

      //calculate the centuries since 2000
      var T;
      if (unix_time){
          mjd = this.getJulianFromUnix(t*1000);
          T = (mjd-2451544.5)/36525;
      }
      else{
          T = (t-2000)/100;
      }

      //Keplerian parameters
      var a = parameters[planet_idx][0]+T*parameters[planet_idx][6];
      var e = parameters[planet_idx][1]+T*parameters[planet_idx][7];
      var I = (parameters[planet_idx][2]+T*parameters[planet_idx][8])*Math.PI/180;
      var l = (parameters[planet_idx][3]+T*parameters[planet_idx][9])*Math.PI/180;
      var w_line = (parameters[planet_idx][4]+T*parameters[planet_idx][10])*Math.PI/180;
      var om = (parameters[planet_idx][5]+T*parameters[planet_idx][11])*Math.PI/180;

      var w = w_line - om;
      var M = l - w_line; //TODO: Add the corrections for jupiter to pluto here.
      M = ( (M + Math.PI) % ( 2*Math.PI) ) - Math.PI;
      var E;

      //Calculate angle
      if (circular===false){
          E = this.keplersEq(M, e);
      }
      else{
          E = M;
      }

      //Get the position in the plane of the ellipse.
      var x = a*(Math.cos(E)-e);
      var y = a*Math.sqrt(1-e*e)*Math.sin(E);

      //Get the position in the ecliptic plane.
      var pos_ecl = [(Math.cos(w)*Math.cos(om)-Math.sin(w)*Math.sin(om)*Math.cos(I))*x + (-Math.sin(w)*Math.cos(om)-Math.cos(w)*Math.sin(om)*Math.cos(I))*y,
          (Math.cos(w)*Math.sin(om)+Math.sin(w)*Math.cos(om)*Math.cos(I))*x + (-Math.sin(w)*Math.sin(om)+Math.cos(w)*Math.cos(om)*Math.cos(I))*y,
          Math.sin(w)*Math.sin(I)*x + Math.cos(w)*Math.sin(I)*y
      ];
      if(icrf === false){
          return pos_ecl;
      }
      else{
          var eps_rad = 23.43928*Math.PI/180;
          var pos_icrf = [pos_ecl[0],
              Math.cos(eps_rad)*pos_ecl[1]-Math.sin(eps_rad)*pos_ecl[2],
              Math.sin(eps_rad)*pos_ecl[1]+Math.cos(eps_rad)*pos_ecl[2]
            ];
          return pos_icrf;
      }
  }
}

/*
//Examples
console.log(planetPosition.getPos(2, 2000, {icrf:true}));
console.log(planetPosition.getPos('jupiter', 2017.101, {circular:true}));
console.log(planetPosition.getPos("earth", 2100.5, {icrf:true, circular:true}));
console.log(planetPosition.getPos(0, 2000.7));
*/
