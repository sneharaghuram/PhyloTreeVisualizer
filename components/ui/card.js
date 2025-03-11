// Card.js
import React from "react";

// Card Component - acts as the container
export const Card = ({ children, className, ...props }) => {
  return (
    <div
      className={`bg-white shadow-md rounded-lg p-4 ${className}`}
      {...props}
    >
      {children}
    </div>
  );
};

// CardHeader - used for the header or title section
export const CardHeader = ({ children, className, ...props }) => {
  return (
    <div className={`border-b-2 pb-4 ${className}`} {...props}>
      {children}
    </div>
  );
};

// CardTitle - used for the title inside the CardHeader
export const CardTitle = ({ children, className, ...props }) => {
  return (
    <h2
      className={`text-xl font-semibold text-gray-800 ${className}`}
      {...props}
    >
      {children}
    </h2>
  );
};

// CardContent - used for the body/content inside the Card
export const CardContent = ({ children, className, ...props }) => {
  return (
    <div className={`mt-4 ${className}`} {...props}>
      {children}
    </div>
  );
};
